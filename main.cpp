#include <fstream>
#include <sstream>
#include <ostream>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <iostream>


#include "knn.h"


typedef CGAL::Gmpq Number_type;
typedef CGAL::Cartesian_d<Number_type> Kernel;
typedef Kernel::Point_d Point_d;

Point_d read_point(size_t d, std::ifstream &is) {

    std::vector<Number_type> point_to_be;
    for (int i = 0; i < d; ++i) {
        Kernel::FT c;
        is >> c;
        point_to_be.push_back(c);
    }
    Point_d point(static_cast<int>(d), point_to_be.begin(), point_to_be.end());
    return point;
}

Kernel::FT squared_dist_d(Point_d p, Point_d q, size_t d) {
    Kernel::FT res = 0;
    for (size_t i = 0; i < d; ++i)
        res += (p[i] - q[i]) * (p[i] - q[i]);

    return res;
}


using namespace std;

//argv: ./a.out dimension input_points_file query_points_file
int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Format: ./a.out dimension input_points_file query_points_file" << std::endl;
        return 1;
    }

    auto d = boost::lexical_cast<size_t>(argv[1]);

    const auto *build_filename = argv[2];
    std::ifstream build_is;
    build_is.open(build_filename);
    if (!build_is.is_open()) {
        std::cerr << "Failed to open " << build_filename << "!" << std::endl;
        return -1;
    }
    size_t bulidn;
    build_is >> bulidn;
    std::vector<Point_d> build_points;

    //cout << "points:" << endl;
    for (size_t i = 0; i < bulidn; ++i) {
        Point_d p  = read_point(d, build_is);
        build_points.push_back(p);
        //cout << i << " " << p << endl;
    }

    const auto *test_filename = argv[3];
    std::ifstream test_is;
    test_is.open(test_filename);
    if (!test_is.is_open()) {
        std::cerr << "Failed to open " << test_filename << "!" << std::endl;
        return -1;
    }
    size_t testn;
    test_is >> testn;
    std::vector<Point_d> test_points;
    //cout << endl << "queries:" << endl;
    for (auto i = 0; i < testn; ++i) {
        Point_d p  = read_point(d, test_is);
        test_points.push_back(p);
        //cout << i << " " << p << endl;
    }
    //cout << endl;
    //create a set of the build points
    std::set<Point_d> pset(build_points.begin(), build_points.end());


    Knn<Kernel> knn(d, build_points.begin(), build_points.end());
    for (size_t k = 2; k < 10; ++k) {
        for (auto it = test_points.begin(); it != test_points.end(); ++it) {
            std::vector<Point_d> res;
            res.reserve(k);
            boost::timer timer;
            knn.find_points(k, *it, std::back_inserter(res));
            auto secs = timer.elapsed();
            std::cout << secs << " time" << std::endl;

           /* cout << "found the following points:" <<endl;
            for(auto &p: res)
                cout << p << endl;
            cout << endl;*/

            //----------------------------------------------//
            //code for testing the output of knn.find_points
            //----------------------------------------------//

            //1. we should get k points or less if k > n
            if ((k != res.size() && k <= build_points.size()) ||
                res.size() > k) {
                std::cerr << "Not enough neighbors" << std::endl;
                return -1;
            }

            //2. make sure that all reported points are points
            //    from the build_points
            for (auto rit = res.begin(); rit != res.end(); ++rit) {
                if (pset.find(*rit) == pset.end()) {
                    std::cerr << "Reported point was not part of the build point set" << std::endl;
                    return -1;
                }
            }

            //3. make sure that each reported point appears once in res
            std::set<Point_d> res_set(res.begin(), res.end());
            if (res_set.size() < res.size()) {
                std::cerr << "A point was reported as a neighbor more than once" << std::endl;
                return -1;

            }


            //4. make sure that there is no point from build_points that is
            //   closer to the query point than the farthest point in res
            //   and that is not reported in res.
            //   Make sure that every point whose distance to the query is as the
            //   distance of the farthest point to the query, is lexicographically greater than the farthest point in res.

            //find maximal in res (point whose dist to test_point it is maximal)

            size_t ind_of_max = 0;
            Kernel::FT maxdist = squared_dist_d(res[0], *it, d);
            for (auto j = 1; j < res.size(); ++j) {
                Kernel::FT j_dist = squared_dist_d(res[j], *it, d);
                if (maxdist < j_dist) {
                    ind_of_max = static_cast<size_t>(j);
                    maxdist = j_dist;
                } else if (maxdist == j_dist) {
                    if (res[ind_of_max] < res[j])
                        ind_of_max = static_cast<size_t>(j);
                }
            }


            Point_d &max_neighbor = res[ind_of_max];

            //verify that no other point in the data set is closer to *it then the farthest in res.
            //If there is a point whose distance from *it is the same then verify it is lexicographically greater than the point in res obtaining the max distance
            for (auto &build_point : build_points) {
                if (build_point == max_neighbor)
                    continue;
                Kernel::FT pit_dist = squared_dist_d(build_point, *it, d);
                if ((pit_dist == maxdist && build_point < max_neighbor) || (pit_dist < maxdist)) {
                    //verify that *pit is reported in res
                    if (res_set.find(build_point) == res_set.end()) {
                        std::cerr << "A potential neighbor was not reported" << build_point << std::endl;
                        return -1;
                    }
                }
            }
            cout << "k " << k << " passed" << endl;

        }
    }

    return 0;
}
