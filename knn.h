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
      /*  for (int i = 1; beginPoints != endPoints; ++beginPoints, ++i)
            this->points.push_back(make_tuple(beginPoints, i, false));
        auto points_index_by_axis = sort_by_each_axis(d);
        this->root = unique_ptr(new KnnNode(points, points_index_by_axis, points_index_by_axis[0], points.size(), 0, d));
        for(int i = 0; i<d; ++i)
            delete[] points_index_by_axis[i];
        delete[] points_index_by_axis;*/
    }


    //input:   const reference to a d dimensional vector which represent a d-point.
    //output:  a vector of the indexes of the k-nearest-neighbors points todo : indexes or point_d??
    template<typename OutputIterator>
    OutputIterator find_points(size_t k, const Point_d &it, OutputIterator oi) {
      /*  KNearestPointRepository knnRepository(points, *root, it, k);
        int i = 0;
        for (int point_index: knnRepository.get_results())
            oi[i++] = points[point_index];*/
        return oi;
    }

private:
  /*  int **sort_by_each_axis(size_t d) {
        int **res = new int *[d];
        for (int i = 0; i < d; ++i) {
            res[i] = new int[points.size()];
            for (int j = 0; j < points.size(); ++j)
                res[i][j] = j;
            sort(res[i], res[i] + points.size(), [](const int &a, const int &b) -> bool {
                return points[a][i] < points[b][i];
            });
        }
        return res;
    }

    class KNearestPointRepository {
        class CompareDistLess{
            bool operator()(pair<int,typename Kernel::FT> &a, pair<Point_d,typename Kernel::FT> &b){
                return a.second<b.second;
            }
        };
        class CompareDistGreater{
            bool operator()(pair<Knn<Kernel>::KnnNode & ,typename  Kernel::FT> &a, pair<Knn<Kernel>::KnnNode &,typename  Kernel::FT> &b){
                return a.second>b.second;
            }
        };

        Kernel::FT squared_dist_d(Point_d p, Point_d q, size_t d)
        {
            Kernel::FT res = 0;
            for (size_t i=0; i < d ;++i)
                res += (p[i] - q[i]) * (p[i] - q[i]);

            return res;
        }

        // todo use another type maybe that can be iterated without pop... cost a lot.
        typedef priority_queue<pair<int,Kerenel::FT>, vector<pair<int,Kerenel::FT>>, CompareDistLess> points_priority_queue_t;
        typedef priority_queue<pair<KnnNode &, Kerene::FT>, vector<pair<KnnNode &, Kerene::FT>>, CompareDistGreater> nodes_priority_queue_t;
    public :
        KNearestPointRepository(vector<Point_d> &points, KnnNode &root, Point_d &current_query_point,
                                size_t k) : points(points), current_query_point(current_query_point),
                                            k(k) {
            points_priority_queue_t knnPointQueue(compare_points);
            nodes_priority_queue_t knnNodeQueue(compare_nodes);
            knnNodeQueue.push(root);
            while (knnNodeQueue.size()) {
                pair<KnnNode, Kernel::FT> currentNode = knnNodeQueue.top();
                knnNodeQueue.pop();
                if (currentNode.first.size() <= k - knnPointQueue.size() ||
                    currentNode.first().size() <= MAX_NUMBER_OF_POINTS_IN_NODE) {
                    for (int point : currentNode.first.get_points()) {
                        Kernel::FT distance = squared_dist_d(points[point], current_query_point);
                        if (knnPointQueue.size() < k)
                            knnPointQueue.push(make_pair(point, distance));
                        else if (knnPointQueue.top().second > distance) {// the points is closer then the father one
                            knnPointQueue.push(make_pair(point,distance));
                            knnPointQueue.pop();
                        }
                    }
                } else {
                    // queue is full.
                    if (currentNode.second > knnPointQueue.top().second) {
                        // if the distance to the most fat away point is smaller then no need to continue;
                        break;
                    } else {
                        closestPoint = currentNode.get_left().get_closest_point_possible(current_query_point);
                        distance = squared_dist_d(closestPoint, current_query_point);
                        nodes_priority_queue.push(make_pair(currentNode.get_left(), distance));
                        closestPoint = currentNode.get_right().get_closest_point_possible(current_query_point);
                        distance = squared_dist_d(closestPoint, current_query_point);
                        nodes_priority_queue.push(make_pair(currentNode.get_right(), distance));
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
    };

    class KnnNode {
    public:
        KnnNode(vector<Point_d> &points, int **points_indexes_sorted_by_axis, const int *points_indexes,
                size_t number_of_points, size_t index_to_sort_by,
                size_t d) : points(points), d(d) {
            if (number_of_points == 0)
                return; // is empty.
            for (int i = 0; i < number_of_points; ++i)
                this->points_indexes.push_back(points_indexes[i]);
            update_min_max(points_indexes_sorted_by_axis, number_of_points);
            if (number_of_points <= MAX_NUMBER_OF_POINTS_IN_NODE)
                return; // no need to split.
            size_t m = number_of_points / 2;

            int **left_points_indexes_ordered_by_axis, **right_points_indexes_ordered_by_axis;

            this->create_splitted_sorted_by_axis_array(points_indexes_sorted_by_axis, &left_points_indexes_ordered_by_axis,
                                                 &right_points_indexes_ordered_by_axis, number_of_points,
                                                 index_to_sort_by);

            this->left =unique_ptr( new KnnNode(points, left_points_indexes_ordered_by_axis,
                               points_indexes_sorted_by_axis[index_to_sort_by], m, (index_to_sort_by + 1) % d, d));
            this->right = unique_ptr(new KnnNode(points, right_points_indexes_ordered_by_axis,
                                points_indexes_sorted_by_axis[index_to_sort_by] + m, number_of_points - m,
                                (index_to_sort_by + 1) % d, d));

            for (int i = 0; i < d; i++)
                delete[] left_points_indexes_ordered_by_axis[i], delete[] right_points_indexes_ordered_by_axis[i];
            delete[] left_points_indexes_ordered_by_axis, delete[] right_points_indexes_ordered_by_axis;
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
        unique_ptr<KnnNode> left, right;

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
            for (int i = 0; i < m; ++i)
                leftIndexesSet.insert(origin[index_to_sort_by][i]);
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

    vector<Point_d> points;
    unique_ptr<KnnNode> root;*/
};