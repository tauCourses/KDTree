#include <vector>
#include <queue>

#include "knnNode.h"
#define MAX_NUMBER_OF_POINTS_IN_NODE 3

using namespace std;

template<typename Kernel>
class Knn {
public:
    typedef typename Kernel::Point_d Point_d;

    // input: dimension - the dimension
    // iterator to the input points
    template<typename InputIterator> Knn(size_t d, InputIterator beginPoints, InputIterator endPoints) :
        dimension(d) {
        for (; beginPoints != endPoints; ++beginPoints)
            this->points.push_back(*beginPoints);
        this->update_sorted_indexes();

        this->root = unique_ptr<KnnNode<Kernel>>(new KnnNode<Kernel>(sorted_indexes, 0, this->dimension));
    }


    //input:   const reference to a dimension dimensional vector which represent a dimension-point.
    //output:  a vector of the indexes of the k-nearest-neighbors points todo : indexes or point_d??
    template<typename OutputIterator>
    OutputIterator find_points(size_t k, const Point_d &it, OutputIterator oi) {
        return this->root->find_points(k, it, oi);
    }

private:
    static bool compareIndexes(const pair<Point_d*, int>& a, const pair<Point_d*, int>& b)
    {
        return (*(a.first))[a.second] < (*(b.first))[b.second];
    }

    void update_sorted_indexes() {
        this->sorted_indexes.resize(this->dimension);
        for (int i = 0; i < dimension; ++i) {
            this->sorted_indexes[i].resize(this->points.size());
            for (int j = 0; j < points.size(); ++j)
                this->sorted_indexes[i][j] = &this->points[j];
            sort(this->sorted_indexes[i].begin(), this->sorted_indexes[i].end(),
                 [&](const Point_d *a, const Point_d *b) -> bool { return (*a)[i] < (*b)[i]; });
/*        for (int i = 0; i < this->dimension; ++i) { TODO->remove, just print
            cout << "for index " << i << ":" << endl;
            for (int j = 0; j < points.size(); ++j)
                cout << get<0>(this->sortedIndexes[i][j]) << endl;
            cout << endl;

        }*/
        }
    }

    size_t dimension;
    vector<Point_d> points;

    vector<vector<Point_d*>> sorted_indexes;
    unique_ptr<KnnNode<Kernel>> root;
};