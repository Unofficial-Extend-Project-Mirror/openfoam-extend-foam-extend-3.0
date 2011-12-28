template<int N, int I> 
class unrollVectorSpaceOps
{
public:

    static const int endLoop = (I < N-1) ? 1 : 0;

    template<class V, class S, class EqOp>
    static inline void eqOpS(V& vs, const S& s, EqOp eo,
                                   unsigned int& offset)
    {
        eo(vs.v_[I+offset], s);
        unrollVectorSpaceOps<endLoop*N, endLoop*(I+1)>
            ::eqOpS(vs, s, eo, offset);
    }

    template<class S, class V, class EqOp>
    static inline void SeqOp(S& s, const V& vs, EqOp eo,
                                   unsigned int& offset)
    {
        eo(s, vs.v_[I+offset]);
        unrollVectorSpaceOps<endLoop*N, endLoop*(I+1)>
            ::SeqOp(s, vs, eo, offset);
    }

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1& vs1, const V2& vs2, EqOp eo,
                                  unsigned int& offset)
    {
        eo(vs1.v_[I+offset], vs2.v_[I+offset]);
        unrollVectorSpaceOps<endLoop*N, endLoop*(I+1)>
            ::eqOp(vs1, vs2, eo, offset);
    }


    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1& vs1, const S& s, Op o,
                                  unsigned int& offset)
    {
        vs.v_[I+offset] = o(vs1.v_[I+offset], s);
        unrollVectorSpaceOps<endLoop*N, endLoop*(I+1)>
            ::opVS(vs, vs1, s, o, offset);
    }

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S& s,  const V1& vs1, Op o,
                                  unsigned int& offset)
    {
        vs.v_[I+offset] = o(s, vs1.v_[I+offset]);
        unrollVectorSpaceOps<endLoop*N, endLoop*(I+1)>
            ::opSV(vs, s, vs1, o, offset);
    }

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1& vs1, const V1& vs2, Op o,
                                unsigned int& offset)
    {
        vs.v_[I+offset] = o(vs1.v_[I+offset], vs2.v_[I+offset]);
        unrollVectorSpaceOps<endLoop*N, endLoop*(I+1)>
            ::op(vs, vs1, vs2, o, offset);
    }

    template<class V, class V1, class V2, class Op>
    static inline void opVV(V& vs, const V1& vs1, const V2& vs2, Op o,
                                  unsigned int& offset)
    {
        vs.v_[I+offset] = o(vs1.v_[I+offset], vs2.v_[I+offset]);
        unrollVectorSpaceOps<endLoop*N, endLoop*(I+1)>
            ::opVV(vs, vs1, vs2, o, offset);
    }

    template<class S, class V1, class V2, class EqOp, class Op>
    static inline void SopEqOpVV(S& s, const V1& vs1, const V2& vs2,
                                       EqOp eo, Op o, unsigned int& offset)
    {
        eo(s, o(vs1.v_[I+offset], vs2.v_[I+offset]));
        unrollVectorSpaceOps<endLoop*N, endLoop*(I+1)>
            ::SopEqOpVV(s, vs1, vs2, eo, o, offset);
    }
};


template<>
class unrollVectorSpaceOps<0, 0>
{
public:

    template<class V, class S, class EqOp>
    static inline void eqOpS(V&, const S&, EqOp,
                                   unsigned int& offset)
    {}

    template<class S, class V, class EqOp>
    static inline void SeqOp(S&, const V&, EqOp,
                                   unsigned int& offset)
    {}

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1&, const V2&, EqOp,
                                  unsigned int& offset)
    {}


    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1&, const S&, Op,
                                  unsigned int& offset)
    {}

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S&, const V1&, Op,
                                  unsigned int& offset)
    {}

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1&, const V1&, Op,
                                unsigned int& offset)
    {}
    
    template<class V, class V1, class V2, class Op>
    static inline void opVV(V& vs, const V1& vs1, const V2& vs2, Op o,
                                  unsigned int& offset)
    {}

    template<class S, class V1, class V2, class EqOp, class Op>
    static inline void SopEqOpVV(S& s, const V1& vs1, const V2& vs2,
                                  EqOp eo, Op o, unsigned int& offset)
    {}
};


template<int N, int I> 
class VectorSpaceOps
{
public:

    static const int unrollN = VECTORSPACE_UNROLLN;
    static const int nOuter = N / unrollN;
    static const int nRemain = N % unrollN;

    template<class V, class S, class EqOp>
    static inline void eqOpS(V& vs, const S& s, EqOp eo)
    {
        unsigned int cmpt=0;
        
        for(int i=0; i<nOuter; i++)
        {
            unrollVectorSpaceOps<unrollN, 0>::eqOpS(vs, s, eo, cmpt);
            cmpt += unrollN;
        }
        
        unrollVectorSpaceOps<nRemain, 0>::eqOpS(vs, s, eo, cmpt);
    }

    template<class S, class V, class EqOp>
    static inline void SeqOp(S& s, const V& vs, EqOp eo)
    {
        unsigned int cmpt=0;
        
        for(int i=0; i<nOuter; i++)
        {
            unrollVectorSpaceOps<unrollN, 0>::SeqOp(s, vs, eo, cmpt);
            cmpt += unrollN;
        }
        
        unrollVectorSpaceOps<nRemain, 0>::SeqOp(s, vs, eo, cmpt);
    }

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1& vs1, const V2& vs2, EqOp eo)
    {
        unsigned int cmpt=0;
        
        for(int i=0; i<nOuter; i++)
        {
            unrollVectorSpaceOps<unrollN, 0>::eqOp(vs1, vs2, eo, cmpt);
            cmpt += unrollN;
        }
        
        unrollVectorSpaceOps<nRemain, 0>::eqOp(vs1, vs2, eo, cmpt);
    }


    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1& vs1, const S& s, Op o)
    {
        unsigned int cmpt=0;
        
        for(int i=0; i<nOuter; i++)
        {
            unrollVectorSpaceOps<unrollN, 0>::opVS(vs, vs1, s, o, cmpt);
            cmpt += unrollN;
        }
        
        unrollVectorSpaceOps<nRemain, 0>::opVS(vs, vs1, s, o, cmpt);
    }

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S& s,  const V1& vs1, Op o)
    {
        unsigned int cmpt=0;
        
        for(int i=0; i<nOuter; i++)
        {
            unrollVectorSpaceOps<unrollN, 0>::opSV(vs, s, vs1, o, cmpt);
            cmpt += unrollN;
        }
        
        unrollVectorSpaceOps<nRemain, 0>::opSV(vs, s, vs1, o, cmpt);
    }

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1& vs1, const V1& vs2, Op o)
    {
        unsigned int cmpt=0;
        
        for(int i=0; i<nOuter; i++)
        {
            unrollVectorSpaceOps<unrollN, 0>::op(vs, vs1, vs2, o, cmpt);
            cmpt += unrollN;
        }
        
        unrollVectorSpaceOps<nRemain, 0>::op(vs, vs1, vs2, o, cmpt);
    }

    template<class V, class V1, class V2, class Op>
    static inline void opVV(V& vs, const V1& vs1, const V2& vs2, Op o)
    {
        unsigned int cmpt=0;
        
        for(int i=0; i<nOuter; i++)
        {
            unrollVectorSpaceOps<unrollN, 0>::opVV(vs, vs1, vs2, o, cmpt);
            cmpt += unrollN;
        }
        
        unrollVectorSpaceOps<nRemain, 0>::opVV(vs, vs1, vs2, o, cmpt);
    }

    template<class S, class V1, class V2, class EqOp, class Op>
    static inline void SopEqOpVV(S& s, const V1& vs1, const V2& vs2, EqOp eo, Op o)
    {
        unsigned int cmpt=0;
        
        for(int i=0; i<nOuter; i++)
        {
            unrollVectorSpaceOps<unrollN, 0>::SopEqOpVV(s, vs1, vs2, eo, o, cmpt);
            cmpt += unrollN;
        }
        
        unrollVectorSpaceOps<nRemain, 0>::SopEqOpVV(s, vs1, vs2, eo, o, cmpt);
    }
};
