#include <iostream>
#include <vector>
#include <cmath>
#include <boost/timer.hpp>

typedef boost::timer timer_type;

typedef double fp_type;
typedef std::vector<fp_type> vector_type;

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
                q[n] = fmod(q[n] + p[n] + pi_2, pi_2);
                p[n] = fmod(p[n] + m_K*sin(q[n]) + pi, pi_2) - pi;
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

    vector_type q(N, pi);
    vector_type p(N);

    for(int n=0; n<N; ++n)
    {
        q[n] = -pi + pi_2*n/N;
    }

    std_map map(0.5);

    timer_type timer;

    map.iterate_N(q, p, steps);

    std::clog << "Iteration finished, " << N << " trajectories, " << steps << " steps: ";
    std::clog << timer.elapsed() << " s" << std::endl;


    for(int n=0; n<N; ++n)
    {
        std::cout << q[n] << '\t' << p[n] << std::endl;
    }
}
