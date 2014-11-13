#include <iostream>
#include <vector>
#include <cmath>
#include <boost/timer.hpp>

#include <boost/simd/sdk/simd/pack.hpp>
#include <boost/simd/sdk/simd/io.hpp>
#include <boost/simd/memory/allocator.hpp>
#include <boost/simd/include/functions/splat.hpp>
#include <boost/simd/include/functions/plus.hpp>
#include <boost/simd/include/functions/multiplies.hpp>
#include <boost/simd/arithmetic/functions.hpp>

#include <nt2/trigonometric/functions.hpp>

namespace simd = boost::simd;

typedef boost::timer timer_type;

typedef double fp_type;
typedef simd::pack<fp_type> simd_pack;
typedef std::vector<simd_pack, simd::allocator<simd_pack> > vector_type;

static const size_t pack_size = simd_pack::static_size;

const fp_type pi = 3.14159265358979312;
const fp_type pi_2 = 2*pi;

struct std_map
{
    const fp_type m_K;

    std_map(const fp_type K) : m_K(K)
    { }

    void iterate_N(vector_type &q, vector_type &p, const int steps)
    {
        const int N = q.size();
        for(int n=0; n<N; ++n)
        {
            for(int step=0; step<steps; step++)
            {
                q[n] = simd::mod(q[n] + p[n] + pi_2, pi_2);
                // put sin in temporary to keep compilation time reasonable...
                simd_pack tmp = m_K*nt2::sin(q[n]);
                p[n] = simd::mod(p[n] + tmp + pi, pi_2) - pi;
            }
        }
    }
};

int main(int argc, char *argv[])
{
    if(argc<3)
    {
        std::cerr << "Expected size and steps as parameter" << std::endl;
        exit(1);
    }
    const size_t N = atoi(argv[1]);
    const size_t steps = atoi(argv[2]);

    vector_type q(N/pack_size);
    vector_type p(N/pack_size);

    // initialize the simd packs
    for(int n=0; n<N/pack_size; ++n)
    {
        for( int j=0 ; j<pack_size; ++j )
            p[n][j] = pi_2*static_cast<fp_type>(n*pack_size+j)/N - pi;
        q[n] = simd::splat<simd_pack>(pi);
    }

    std_map map(0.5);

    timer_type timer;

    map.iterate_N(q, p, steps);

    std::clog << "Iteration finished, " << N << " trajectories, " << steps << " steps: ";
    std::clog << timer.elapsed() << " s" << std::endl;

    // print results
    for(int n=0; n<N/pack_size; ++n)
    {
        for( int j=0 ; j<pack_size; ++j )
            std::cout << q[n][j] << '\t' << p[n][j] << std::endl;
    }
}
