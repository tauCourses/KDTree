#include <vector>
#include <queue>

#define MAX_NUMBER_OF_POINTS_IN_NODE 3

using namespace std;

template<typename Kernel>
class Knn {
public:
    typedef typename Kernel::Point_d Point_d;

    // input: d - the dimension
    // iterator to the input points
    template<typename InputIterator>
    Knn(size_t d, InputIterator beginPoints, InputIterator endPoints) {
        for (int i = 1; beginPoints != endPoints; ++beginPoints, ++i)
            points.push_back(make_tuple(beginPoints, i, false));
        auto points_index_by_axis = sort_by_each_axis(d);
        root = unique_ptr(new KnnNode(points, points_index_by_axis, points_index_by_axis[0], points.size(), 0, d));
    }


    //input:   const reference to a d dimensional vector which represent a d-point.
    //output:  a vector of the indexes of the k-nearest-neighbors points todo : indexes or point_d??
    template<typename OutputIterator>
    OutputIterator find_points(size_t k, const Point_d &it, OutputIterator oi) {
        KNearestPointRepository knnRepository(points, *root, it, k);
        for (int point_index: knnRepository.get_results()) {
            //if need to return pointers:
            oi++ = point_index; // todo check if indexes are 0 based
            //else
            oi++ = points[point_index];
        }
        return oi;
    }

private:
    vector<Point_d> points;
    unique_ptr<KnnNode> root;

    int **sort_by_each_axis(size_t d) {
        int **res = new int *[d];
        for (int i = 0; i < d; ++i) {
            res[i] = new int[points.size()];
            for (int j = 0; j < points.size(); ++j) res[i][j] = j;
            sort(res[i], res[i] + points.size(), [](const int &a, const int &b) -> bool {
                return points[a][i] < points[b][i];
            });
        }
        return res;
    }

    class KNearestPointRepository {
        // todo use another type maybe that can be iterated without pop... cost a lot.
        typedef priority_queue<int, vector<int>, function<bool(int, int)>> points_priority_queue;
        typedef priority_queue<Knn<Kernel>::KnnNode &, vector<Knn<Kernel>::KnnNode &>, std::function<bool(
                Knn<Kernel>::KnnNode &, Knn<Kernel>::KnnNode &)>> nodes_priority_queue;
    public :
        KNearestPointRepository(vector<Point_d> &points, KnnNode &root, Point_d &current_query_point,
                                size_t k) : points(points), current_query_point(current_query_point),
                                            k(k) {
            points_priority_queue knnPointQueue(compare_points);
            nodes_priority_queue knnNodeQueue(compare_nodes);
            knnNodeQueue.push(root);
            while (knnNodeQueue.size()) {
                KnnNode &currentNode = knnNodeQueue.top();
                knnNodeQueue.pop();
                if (currentNode.size() <= k - knnPointQueue.size() ||
                    currentNode.size() <= MAX_NUMBER_OF_POINTS_IN_NODE) {
                    for (int point : currentNode.get_points())
                        knnPointQueue.push(point);
                } else {
                    //queue is full.
                    Point_d closestPoint = currentNode.get_closest_point_possible(current_query_point);
                    if (0) { // if the distance to the most fat away point is smaller then no need to continue;
                        break;
                    } else {
                        nodes_priority_queue.push(currentNode.get_left());
                        nodes_priority_queue.push(currentNode.get_right());
                    }
                }
            }
            while (knnPointQueue.size()) {
                result.push_back(knnPointQueue.top());
                knnPointQueue.pop();
            }
        }

        const vector<int> &get_results() const { return result; }

    private:
        const vector<const Point_d> &points;
        const Point_d &current_query_point;
        size_t k;
        vector<int> result;

        bool compare_points(const int lhs_point_index, const int rhs__point_index) const {
            // todo now return which of the points is the closest one to current_query_point.
            return 0;
        }

        bool compare_nodes(const KnnNode &lhs, const KnnNode &rhs) const {
            Point_d closest_point_in_lhs = lhs.get_closest_point_possible(current_query_point);
            Point_d closest_point_in_rhs = rhs.get_closest_point_possible(current_query_point);
            // todo now return which of the points is the MOST FAR AWAY one to current_query_point.
            return 0;
        }
    };

    class KnnNode {
    public:
        KnnNode(vector<Point_d> &points, int **points_indexes_sorted_by_axis, const int *points_indexes,
                size_t number_of_points, size_t index_to_sort_by,
                size_t d) : points(points), d(d) {
            if (number_of_points == 0) return; // is empty.
            for (int i = 0; i < number_of_points; ++i)
                this->points_indexes.push_back(points_indexes[i]);
            update_min_max(points_indexes_sorted_by_axis, number_of_points);
            if (number_of_points <= MAX_NUMBER_OF_POINTS_IN_NODE) return; // no need to split.
            size_t m = number_of_points / 2;

            int **left_points_indexes_ordered_by_axis, **right_points_indexes_ordered_by_axis;

            create_splitted_sorted_by_axis_array(points_indexes_sorted_by_axis, &left_points_indexes_ordered_by_axis,
                                                 &right_points_indexes_ordered_by_axis, number_of_points,
                                                 index_to_sort_by);

            left = new KnnNode(points, left_points_indexes_ordered_by_axis,
                               points_indexes_sorted_by_axis[index_to_sort_by], m, (index_to_sort_by + 1) % d, d);
            right = new KnnNode(points, right_points_indexes_ordered_by_axis,
                                points_indexes_sorted_by_axis[index_to_sort_by] + m, number_of_points - m,
                                (index_to_sort_by + 1) % d, d);

            for (int i = 0; i < d; i++)
                delete[] left_points_indexes_ordered_by_axis[i], delete[] right_points_indexes_ordered_by_axis[i];
            delete[] left_points_indexes_ordered_by_axis, delete[] right_points_indexes_ordered_by_axis;
        }

        ~KnnNode() {
            delete left;
            delete right;
        }

        // todo make this faster by caching the last result.
        Point_d get_closest_point_possible(Point_d &point) const {
            Point_d res;
            for (int i = 0; i < d; i++) {
                if (min_point[i] > point[i]) {
                    res[i] = min_point[i];
                } else if (max_point[i] < point[i]) {
                    res[i] = max_point[i];
                } else {
                    res[i] = point[i];
                }
            }
            return res;
        }

        size_t size() const {
            return points_indexes.size();
        }

        const vector<int> &get_points() const {
            return points_indexes;
        }

        KnnNode &get_left() const {
            return *left;
        }

        KnnNode &get_right() const {
            return *right;
        }

    private:
        vector<Point_d> &points;
        vector<int> points_indexes;
        size_t d;
        Point_d min_point;
        Point_d max_point;
        KnnNode *left = nullptr, *right = nullptr;

        void update_min_max(int **points_indexes_ordered_by_axis, size_t number_of_points) {
            for (size_t i = 0; i < d; i++) {
                min_point[i] = points[points_indexes_ordered_by_axis[i][0]][i];
                max_point[i] = points[points_indexes_ordered_by_axis[i][number_of_points - 1]][i];
            }
        }

        void create_splitted_sorted_by_axis_array(int **origin, int ***pleft, int ***pright, size_t n,
                                                  size_t index_to_sort_by) {
            size_t m = n / 2;
            set<int> leftIndexesSet;
            for (int i = 0; i < m; ++i) leftIndexesSet.insert(origin[index_to_sort_by][i]);
            (*pleft) = new int *[d], (*pright) = new int *[d];
            for (int i = 0; i < d; ++i) {
                (*pleft)[i] = new int[m], (*pright)[i] = new int[n - m];
                int left_index = 0, right_index = 0;
                for (int j = 0; j < m; ++j) {
                    if (leftIndexesSet.find(origin[i][j]) != leftIndexesSet.end())
                        (*pleft)[i][left_index++] = origin[i][j];
                    else
                        (*pright)[i][right_index++] = origin[i][j];
                }
            }
        }
    };
};