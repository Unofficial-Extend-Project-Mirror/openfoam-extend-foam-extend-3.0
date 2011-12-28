template<int N, int I> 
class VectorSpaceOps
{
public:

    template<class V, class S, class EqOp>
    static inline void eqOpS(V& vs, const S& s, EqOp eo)
    {
        for(register int i=I; i<N; i++)
        {
            eo(vs.v_[i], s);
        }
    }

    template<class S, class V, class EqOp>
    static inline void SeqOp(S& s, const V& vs, EqOp eo)
    {
        for(register int i=I; i<N; i++)
        {
            eo(s, vs.v_[i]);
        }
    }

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1& vs1, const V2& vs2, EqOp eo)
    {
        for(register int i=I; i<N; i++)
        {
            eo(vs1.v_[i], vs2.v_[i]);
        }
    }

    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1& vs1, const S& s, Op o)
    {
        for(register int i=I; i<N; i++)
        {
            vs.v_[i] = o(vs1.v_[i], s);
        }
    }

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S& s,  const V1& vs1, Op o)
    {
        for(register int i=I; i<N; i++)
        {
            vs.v_[i] = o(s, vs1.v_[i]);
        }
    }

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1& vs1, const V1& vs2, Op o)
    {
        for(register int i=I; i<N; i++)
        {
            vs.v_[i] = o(vs1.v_[i], vs2.v_[i]);
        }
    }

    template<class V, class V1, class V2, class Op>
    static inline void opVV(V& vs, const V1& vs1, const V2& vs2, Op o)
    {
        for(register int i=I; i<N; i++)
        {
            vs.v_[i] = o(vs1.v_[i], vs2.v_[i]);
        }
    }

    template<class S, class V1, class V2, class EqOp, class Op>
    static inline void SopEqOpVV(S& s, const V1& vs1, const V2& vs2, EqOp eo, Op o)
    {
        for(register int i=I; i<N; i++)
        {
            eo(s, o(vs1.v_[i], vs2.v_[i]));
        }
    }
};

template<> 
class VectorSpaceOps<0, 0>
{
public:

    template<class V, class S, class EqOp>
    static inline void eqOpS(V& vs, const S& s, EqOp eo)
    {}

    template<class S, class V, class EqOp>
    static inline void SeqOp(S& s, const V& vs, EqOp eo)
    {}

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1& vs1, const V2& vs2, EqOp eo)
    {}

    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1& vs1, const S& s, Op o)
    {}

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S& s,  const V1& vs1, Op o)
    {}

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1& vs1, const V1& vs2, Op o)
    {}

    template<class V, class V1, class V2, class Op>
    static inline void opVV(V& vs, const V1& vs1, const V2& vs2, Op o)
    {}

    template<class S, class V1, class V2, class EqOp, class Op>
    static inline void SopEqOpVV(S& s, const V1& vs1, const V2& vs2, EqOp eo, Op o)
    {}
};