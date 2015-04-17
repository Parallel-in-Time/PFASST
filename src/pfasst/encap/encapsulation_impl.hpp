#include "pfasst/encap/encapsulation.hpp"

#include "pfasst/globals.hpp"


namespace pfasst
{
  namespace encap
  {
    template<typename time>
    Encapsulation<time>::~Encapsulation()
    {}

    /**
     * @throws NotImplementedYet This function is required by Encapsulation
     */
    template<typename time>
    void Encapsulation<time>::zero()
    {
      throw NotImplementedYet("encap");
    }

    /**
     * @throws NotImplementedYet This function is required by Encapsulation
     */
    template<typename time>
    void Encapsulation<time>::copy(shared_ptr<const Encapsulation<time>>)
    {
      throw NotImplementedYet("encap");
    }

    /**
     * @throws NotImplementedYet This function is required by Encapsulation
     */
    template<typename time>
    time Encapsulation<time>::norm0() const
    {
      throw NotImplementedYet("encap");
    }

    /**
     * @throws NotImplementedYet This function is required by Encapsulation
     */
    template<typename time>
    void Encapsulation<time>::saxpy(time a, shared_ptr<const Encapsulation<time>> x)
    {
      UNUSED(a); UNUSED(x);
      throw NotImplementedYet("encap");
    }

    /**
     * @internals
     * The matrix-vector multiplication is implemented using Encapsulation::saxpy() for all non-zero
     * entries of @p mat.
     * @endinternals
     *
     * @todo Consider asserting size of matrix matches sizes of vectors.
     */
    template<typename time>
    void Encapsulation<time>::mat_apply(vector<shared_ptr<Encapsulation<time>>> dst,
                                        time a, Matrix<time> mat,
                                        vector<shared_ptr<Encapsulation<time>>> src,
                                        bool zero)
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      if (zero) {
        for (auto elem : dst) { elem->zero(); }
      }

      for (size_t n = 0; n < ndst; n++) {
        for (size_t m = 0; m < nsrc; m++) {
          auto s = mat(n, m);
          if (s != 0.0) {
            dst[n]->saxpy(a*s, src[m]);
          }
        }
      }
    }

    template<typename time>
    void Encapsulation<time>::post(ICommunicator* comm, int tag)
    {
      UNUSED(comm); UNUSED(tag);
    }

    /**
     * @throws NotImplementedYet This function is required by PFASST
     */
    template<typename time>
    void Encapsulation<time>::send(ICommunicator* comm, int tag, bool blocking)
    {
      UNUSED(comm); UNUSED(tag); UNUSED(blocking);
      throw NotImplementedYet("pfasst");
    }

    /**
     * @throws NotImplementedYet This function is required by PFASST
     */
    template<typename time>
    void Encapsulation<time>::recv(ICommunicator* comm, int tag, bool blocking)
    {
      UNUSED(comm); UNUSED(tag); UNUSED(blocking);
      throw NotImplementedYet("pfasst");
    }

    /**
     * @throws NotImplementedYet This function is required by PFASST
     */
    template<typename time>
    void Encapsulation<time>::broadcast(ICommunicator* comm)
    {
      UNUSED(comm);
      throw NotImplementedYet("pfasst");
    }
  }  // ::pfasst::encap
}  // ::pfasst
