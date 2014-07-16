/*
 * Sweeper for scalar test equation
 */
 
 #ifndef _SCALAR_SWEEPER_HPP_
 #define _SCALAR_SWEEPER_HPP_
 
 #include <complex>
 #include <vector>
 #include <pfasst/encap/imex_sweeper.hpp>
 
 using namespace std;
 
 template<typename time = pfasst::time_precision>
 class ScalarSweeper : public pfasst::encap::IMEXSweeper<time>
 {
 
    typedef pfasst::encap::Encapsulation<time> Encapsulation;
    typedef pfasst::encap::VectorEncapsulation<complex<double>> DVectorT;
    
    public:
        ScalarSweeper(complex<double> lambda, complex<double> y0)
        {
            this->_lambda  = lambda;
            this->_y0      = y0;
            this->_nf1eval = 0;
            this->_nf2eval = 0;
            this->_nf2comp = 0;
        }
        
        ~ScalarSweeper()
        {
            cout << "Number of calls to f1eval: " << this->_nf1eval << endl;
            cout << "Number of calls to f2eval: " << this->_nf2eval << endl;
            cout << "Number of calls to f2comp: " << this->_nf2comp << endl;
        }
        
        void echo_error(time t, bool predict = false)
        {
            int lastnode = this->get_nodes().size()-1;
            shared_ptr<DVectorT> qend = dynamic_pointer_cast<DVectorT>(this->get_state(lastnode));
            assert(qend);
            
            shared_ptr<DVectorT> qex = make_shared<DVectorT>(qend->size());
                       
            this->exact(qex, t);
            double max_err = abs(qend->data()[0] - qex->data()[0])/abs(qex->data()[0]);
            cout << "err: " << scientific << max_err << endl;
        }
        
        void predict(time t, time dt, bool initial)
        {
            pfasst::encap::IMEXSweeper<time>::predict(t, dt, initial);
            this->echo_error(t + dt, true);
        }
        
        void sweep(time t, time dt)
        {
            pfasst::encap::IMEXSweeper<time>::sweep(t, dt);
            this->echo_error(t + dt);
        }
        
        // Here come the functions that implement the imex_sweeper (?) interface; in
        // a first step, the pointer to a general encapsulation object are cast on
        // pointers to the actually used solution structure of type DVectorT
        
        void exact(shared_ptr<Encapsulation> q, time t)
        {
            shared_ptr<DVectorT> q_cast = dynamic_pointer_cast<DVectorT>(q);
            assert(q_cast);
            this->exact(q_cast, t);
        }
                
        void f1eval(shared_ptr<Encapsulation> f, shared_ptr<Encapsulation> q, time t)
        {
            shared_ptr<DVectorT> f_cast = dynamic_pointer_cast<DVectorT>(f);
            assert(f);
            shared_ptr<DVectorT> q_cast = dynamic_pointer_cast<DVectorT>(q);
            assert(q_cast);
            
            this->f1eval(f_cast, q_cast, t);
        }
        
        void f2eval(shared_ptr<Encapsulation> f, shared_ptr<Encapsulation> q, time t)
        {
            shared_ptr<DVectorT> f_cast = dynamic_pointer_cast<DVectorT>(f);
            assert(f);
            shared_ptr<DVectorT> q_cast = dynamic_pointer_cast<DVectorT>(q);
            assert(q_cast);
            
            this->f2eval(f_cast, q_cast, t);
        }
     
        
        void f2comp(shared_ptr<Encapsulation> f, shared_ptr<Encapsulation> q, time t, time dt,
                    shared_ptr<Encapsulation> rhs)
        {
              shared_ptr<DVectorT> f_cast   = dynamic_pointer_cast<DVectorT>(f);
              assert(f_cast);
              shared_ptr<DVectorT> q_cast   = dynamic_pointer_cast<DVectorT>(q);
              assert(q_cast);
              shared_ptr<DVectorT> rhs_cast = dynamic_pointer_cast<DVectorT>(rhs);
              assert(rhs_cast);

              this->f2comp(f_cast, q_cast, t, dt, rhs_cast);            
        }
        
        // Now come the actual implementations working on the DVectorT solutions
           
        void exact(shared_ptr<DVectorT> q, time t)
        {
            q->data()[0] = _y0*exp(_lambda*t);
        }
              
        void f1eval(shared_ptr<DVectorT> f, shared_ptr<DVectorT> q, time t)
        {
        
            // f1 = multiply with imaginary part of lambda
            f->data()[0] = i_complex*imag(this->_lambda)*q->data()[0];
            this->_nf1eval++;
        }
        

        
        void f2eval(shared_ptr<DVectorT> f, shared_ptr<DVectorT> q, time t)
        {
        
            // f2 = multiply with real part of lambda
            f->data()[0] = real(this->_lambda)*q->data()[0];
            this->_nf2eval++;
        }        

        
        void f2comp(shared_ptr<DVectorT> f, shared_ptr<DVectorT> q, time t, time dt,
                    shared_ptr<DVectorT> rhs)
        {
        
            // invert f2=multiply with inverse of real part of lambda
            double inv = 1/(1 - double(dt)*real(this->_lambda));
            //cout << q->data()[0] << endl;
            //cout << inv << endl;
            q->data()[0] = inv*rhs->data()[0];
            f->data()[0] = real(this->_lambda)*q->data()[0];
            this->_nf2comp++;
        }
        
    private:
    
        complex<double> _lambda, _y0;
        int _nf1eval, _nf2eval, _nf2comp;
        const complex<double> i_complex = complex<double>(0,1);
         
 };
 #endif