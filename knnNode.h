#ifndef KDTREE_KNNNODE_H
#define KDTREE_KNNNODE_H

#include <vector>
#include <memory>

using namespace std;

#define MAX_NODE_SIZE 3
template<typename Kernel>
class KnnNode {
public:
    typedef typename Kernel::Point_d Point_d;
    typedef typename Kernel::FT FT;

    explicit KnnNode(vector<vector<Point_d*>>& sorted_indexes, int currentInx, size_t d) :
            sortedIndexes(sorted_indexes), currentInx(currentInx), dimension(d), numberOfPoints(sorted_indexes[0].size())
    {
        if(this->numberOfPoints <= MAX_NODE_SIZE)
        {
            this->atomic = true;
        } else{
            this->initLeftAndRightSortedIndexes();
            this->left = unique_ptr<KnnNode<Kernel>>(new KnnNode<Kernel>(leftSortedIndexes,
                                                                         static_cast<int>((currentInx + 1) % d),
                                                                         this->dimension));
            this->right = unique_ptr<KnnNode<Kernel>>(new KnnNode<Kernel>(rightSortedIndexes,
                                                                         static_cast<int>((currentInx + 1) % d),
                                                                          this->dimension));
        }

    }

    FT minDistance(const Point_d &p)
    {
        FT result = 0;
        for (int i = 0; i < this->dimension; i++) {
            if ((*sortedIndexes[i][numberOfPoints-1])[i] < p[i]) {
                result += ((*sortedIndexes[i][numberOfPoints - 1])[i] - p[i]) *
                          ((*sortedIndexes[i][numberOfPoints - 1])[i] - p[i]);
            }else if ((*sortedIndexes[i][0])[i] > p[i]) {
                result += ((*sortedIndexes[i][0])[i] - p[i]) *
                          ((*sortedIndexes[i][0])[i] - p[i]);
            }
        }
        return result;
    }

    template<typename OutputIterator>
    OutputIterator find_points(size_t k, const Point_d &it, OutputIterator oi) {
        ksQueue best_ks;
        this->getMinPoints(k,it,best_ks);
        vector<Point_d> temp;
        temp.resize(k);
        for(auto i= static_cast<int>(k - 1); i>=0; i--, oi++)
        {
            oi = *(best_ks.top().first);
            best_ks.pop();
        }
        return oi;
    }


private:
    int currentInx;
    bool atomic = false;
    size_t dimension;
    int numberOfPoints;
    vector<vector<Point_d*>>& sortedIndexes;
    vector<vector<Point_d*>> leftSortedIndexes,rightSortedIndexes;
    unique_ptr<KnnNode<Kernel>> left, right;

    class distanceComperator{
    public:
        bool operator()(const pair<Point_d*,FT>& a, const pair<Point_d*,FT>& b) {
            if(a.second > b.second)
                return false;
            else if(a.second == b.second)
                return *a.first < *b.first;
            return true;
        }
    };

    typedef priority_queue<pair<Point_d*, FT>, vector<pair<Point_d*, FT>>,
            typename KnnNode<Kernel>::distanceComperator> ksQueue;

    FT pointDistance(const Point_d& p, Point_d* q) {
        FT res = 0;
        for (size_t i = 0; i < this->dimension; ++i)
            res += (p[i] - (*q)[i]) * (p[i] - (*q)[i]);

        return res;
    }

    void getMinPoints(size_t k, const Point_d &it, ksQueue& best_ks)
    {
        if(atomic)
        {
            for(auto &t:this->sortedIndexes[0])// point in points
            {
                auto distance = pointDistance(it, t);
                if(best_ks.size() < k || distance <= best_ks.top().second) {
                    best_ks.push(pair<Point_d*, FT>(t, distance));
                    if(best_ks.size() > k)
                        best_ks.pop();
                }
            }
        } else{
            FT leftDistance = this->left->minDistance(it);
            FT rightDistance = this->right->minDistance(it);

            if(leftDistance < rightDistance) {
                if (best_ks.size() < k || leftDistance <= best_ks.top().second) {
                    this->left->getMinPoints(k, it, best_ks);
                    if (best_ks.size() < k || rightDistance <= best_ks.top().second)
                        this->right->getMinPoints(k, it, best_ks);
                }
            } else if (best_ks.size() < k || rightDistance <= best_ks.top().second){
                this->right->getMinPoints(k,it,best_ks);
                if(best_ks.size() < k || rightDistance <= best_ks.top().second)
                    this->left->getMinPoints(k, it, best_ks);
            }
        }
    }
    void initLeftAndRightSortedIndexes()
    {
        int medianIndex = this->numberOfPoints / 2;
        set<Point_d*> leftPart;

        for (int i = 0; i < medianIndex; i++) {
            leftPart.insert(this->sortedIndexes[this->currentInx][i]);
        }

        this->leftSortedIndexes.resize(this->dimension);
        this->rightSortedIndexes.resize(this->dimension);
        for (int i = 0; i < this->dimension; i++) {
            for (int j = 0; j < this->numberOfPoints; j++) {
                if (leftPart.find(this->sortedIndexes[i][j]) != leftPart.end())
                    this->leftSortedIndexes[i].push_back(this->sortedIndexes[i][j]);
                else
                    this->rightSortedIndexes[i].push_back(this->sortedIndexes[i][j]);


            }
        }
    }
};



#endif //KDTREE_KNNNODE_H
