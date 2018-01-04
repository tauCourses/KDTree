#include <vector>

template<typename Kernel> class Knn {
public:
    typedef typename Kernel::Point_d Point_d;

    //input: d - the dimension
    // itreator to the input points
    template<typename InputIterator> Knn(size_t d, InputIterator beginPoints, InputIterator endPoints) {
        //your build code
    }

    //input:   const reference to a d dimensional vector which represent a d-point.
    //output:  a vector of the indexes of the k-nearest-neighbors points
    template<typename OutputIterator> OutputIterator find_points(size_t k, const Point_d &it, OutputIterator oi) {
        //your query code
        return oi;
    }
};
