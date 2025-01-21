/*
  Copyright 2015 IRIS AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef DUNE_COMMUNICATOR_HEADER_INCLUDED
#define DUNE_COMMUNICATOR_HEADER_INCLUDED

#include <cassert>
#include <algorithm>
#include <vector>
#include <set>
#include <map>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/communication.hh>

// the following implementation is only available in case MPI is available
#if HAVE_MPI
#include <dune/common/parallel/mpicommunication.hh>
#endif

#include <opm/grid/grid_size.h>

namespace Dune
{
  class SimpleMessageBuffer
  {
    typedef std::vector< char >  BufferType;

    mutable BufferType buffer_;
    const double factor_;
    mutable size_t pos_;
public:
    /** \brief constructor taking memory reserve estimation factor (default is 1.1, i.e. 10% over estimation )
     */
    SimpleMessageBuffer( const double factor = 1.1 )
      : buffer_(), factor_( factor )
    {
      resetReadPosition();
    }

    /** \brief clear the buffer */
    void clear() { buffer_.clear(); resetReadPosition(); }
    /** \brief reset read position of buffer to beginning */
    void resetReadPosition() { pos_ = 0 ; }
    /** \brief return size of buffer */
    size_t size() const { return buffer_.size(); }

    /** \brief reserve memory for 'size' entries  */
    void reserve( const size_t size )
    {
        buffer_.reserve( size );
    }

    /** \brief resize buffer to 'size' entries  */
    void resize( const size_t size )
    {
        buffer_.resize( size );
    }

    /** \brief write value to buffer, value must implement the operator= correctly (i.e. no internal pointers etc.) */
    template <class T>
    void write( const T& value )
    {
      // union to access bytes in value
      const size_t tsize = sizeof( T );
      size_t pos  = buffer_.size();
      const size_t sizeNeeded = pos + tsize ;
      // reserve with some 10% overestimation
      if( buffer_.capacity() < sizeNeeded )
      {
        reserve( size_t(factor_ * sizeNeeded) ) ;
      }
      // resize to size need to store value
      buffer_.resize( sizeNeeded );
      // copy value to buffer
      std::copy_n( reinterpret_cast<const char *> (&value), tsize, buffer_.data()+pos );
    }

    void write( const std::string& str)
    {
        int size = str.size();
        write(size);
        for (int k = 0; k < size; ++k) {
            write(str[k]);
        }
    }

    /** \brief read value from buffer, value must implement the operator= correctly (i.e. no internal pointers etc.) */
    template <class T>
    void read( T& value ) const
    {
      // read bytes from stream and store in value
      const size_t tsize = sizeof( T );
      assert( pos_ + tsize <= buffer_.size() );
      std::copy_n( buffer_.data()+pos_, tsize, reinterpret_cast<char *> (&value) );
      pos_ += tsize;
    }

    void read( std::string& str) const
    {
        int size = 0;
        read(size);
        str.resize(size);
        for (int k = 0; k < size; ++k) {
            read(str[k]);
        }
    }

    /** \brief return pointer to buffer and size for use with MPI functions */
    std::pair< char* , int > buffer() const
    {
      return std::make_pair( buffer_.data(), int(buffer_.size()) );
    }
  };

  /** \brief Point-2-Point communicator for exchange messages between processes */
  template < class MsgBuffer >
  class Point2PointCommunicator : public CollectiveCommunication< MPIHelper::MPICommunicator >
  {
  public:
    /** \brief type of MPI communicator, either MPI_Comm or NoComm as defined in MPIHelper */
    typedef MPIHelper::MPICommunicator  MPICommunicator ;

    /** \brief type of message buffer used */
    typedef MsgBuffer MessageBufferType ;

  protected:
    using BaseType = Dune::Communication<MPICommunicator>;
    typedef Point2PointCommunicator< MessageBufferType > ThisType;

    // starting message tag
    static const int messagetag = 234;

    typedef std::map< int, int > linkage_t;
    typedef std::vector< int >   vector_t;

    linkage_t  sendLinkage_ ;
    linkage_t  recvLinkage_ ;

    vector_t   sendDest_ ;
    vector_t   recvSource_ ;

    mutable vector_t   _recvBufferSizes;
    mutable bool       _recvBufferSizesComputed;

  public :
    using BaseType :: rank;
    using BaseType :: size;

    /* \brief data handle interface that needs to be implemented for use with some of
     * the exchange methods */
    class DataHandleInterface
    {
    protected:
      DataHandleInterface () {}
    public:
      virtual ~DataHandleInterface () {}
      virtual void   pack( const int link, MessageBufferType& os ) = 0 ;
      virtual void unpack( const int link, MessageBufferType& os ) = 0 ;
      // should contain work that could be done between send and receive
      virtual void localComputation () {}
    };

  public:
    /** \brief constructor taking mpi communicator */
    Point2PointCommunicator( const MPICommunicator& mpiComm = MPIHelper::getCommunicator() )
      : BaseType( mpiComm ) { removeLinkage(); }

    /** \brief constructor taking collective communication */
    Point2PointCommunicator( const BaseType& comm ) : BaseType( comm ) { removeLinkage(); }


    /** \brief insert communication request with a set os ranks to send to and a set of ranks to receive from */
    inline void insertRequest( const std::set< int >& sendLinks, const std::set< int >& recvLinks );

    /** \brief return number of processes we will send data to */
    inline int sendLinks () const { return sendLinkage_.size(); }

    /** \brief return number of processes we will receive data from */
    inline int recvLinks () const { return recvLinkage_.size(); }

    /** \brief return vector containing possible recv buffer sizes */
    const vector_t& recvBufferSizes() const { return _recvBufferSizes; }

    /** \brief return send link number for a given send rank number */
    inline int sendLink (const int rank) const
    {
      assert (sendLinkage_.end () != sendLinkage_.find (rank)) ;
      return (* sendLinkage_.find (rank)).second ;
    }

    /** \brief return recv link number for a given recv rank number */
    inline int recvLink (const int rank) const
    {
      assert (recvLinkage_.end () != recvLinkage_.find (rank)) ;
      return (* recvLinkage_.find (rank)).second ;
    }

    /** \brief return vector containing all process numbers we will send to */
    const std::vector< int > &sendDest   () const { return sendDest_; }
    /** \brief return vector containing all process numbers we will receive from */
    const std::vector< int > &recvSource () const { return recvSource_; }

    /** \brief remove stored linkage */
    inline void removeLinkage () ;

    /** \brief exchange message buffers with peers defined by inserted linkage */
    virtual std::vector< MessageBufferType > exchange (const std::vector< MessageBufferType > &) const;

    /** \brief exchange data with peers, handle defines pack and unpack of data */
    virtual void exchange ( DataHandleInterface& ) const;

    /** \brief exchange data with peers, handle defines pack and unpack of data,
     *  if receive buffers are known from previous run and have not changed
     *  communication could be faster */
    virtual void exchangeCached ( DataHandleInterface& ) const;

  protected:
    inline void computeDestinations( const linkage_t& linkage, vector_t& dest );

    // return new tag number for the exchange messages
    int getMessageTag(const grid_size_t increment) const
    {
      const int retTag = this->tag_;
      this->generateNextMessageTag(increment);
      return retTag;
    }

    // return new tag number for the exchange messages
    int getMessageTag() const
    {
      return this->getMessageTag(1u);
    }

  private:
    mutable int tag_{messagetag + 2};

    void generateNextMessageTag(const grid_size_t increment) const
    {
      this->tag_ += increment;

      // Reset to initial value if next tag exceeds MPI standard's maximum
      // message tag guarantee of 2^15-1.
      if (this->tag_ >= (1 << 15) - 1)
      {
        this->tag_ = messagetag + 2;
      }
    }
  };

} // namespace Dune

// include inline implementation
#include "p2pcommunicator_impl.hh"

#endif // #ifndef DUNE_COMMUNICATOR_HEADER_INCLUDED
