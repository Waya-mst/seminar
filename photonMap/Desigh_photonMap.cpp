#include <いろいろ>

define もろもろ

// 返すのは0.0~1.0までの値
inline double rand01() {return (double)rand()/RAND_MAX;}

struct Vec{
    double x{};
    double y{};
    double z{};

    //コンストラクタとかオペレータとか
};

typedef Vec Color;
const Color BackGroundColor(0.0, 0.0, 0.0);

struct Ray{
    // |dir| = 1
    Vec org, dir;
    Ray(const Vec org_, const Vec dir_):
    org(org_), dir(dir_){}
};

enum ReflectionType{
    DIFFUSE,
    SPECULAR,
    REFRACTION,
};

struct Sphere{
    double radius;
    Vec pos;
    Color emission, color;
    ReflectionType ref_type;

    // コンストラクタ

    const double intersect(const Ray &ray){
        // 交差判定
    }
};

static std::vector<Sphere*> spheres = {
    // シーンデータ
}
const int LightID = 0;

struct Photon{
    Vec pos;
    Color power;
    Vec incident;

    Photon(const Vec& pos_, const Color& power_, const Vec& incident_)
    : pos(pos_), power(power_), incident(incident_) {}
};

template<typename T>
class KDTree{
public:
    struct Query{　// 近傍探索する時の条件
        Vec search_position; // 探索中心
        Vec normal; // 探索中心における法線
        double max_distance2; // 探索の最大半径
        size_t max_search_num; // 最大探索点数
    };

    struct ElementForQueue{
        // 
    };

    typedef std::priority_queue<ElementForQueue, std::vector<ElementForQueue>> ResultQueue;

private:
    std::vector<T> points;
    struct KDTreeNode{
        //
    };

    KDTreeNode* root;
    
    // デストラクタ

    void locate_points(){
        // KDTreeの中からクエリの条件に合うフォトンを探索してResultQueueに詰める
        // 実質KD-Treeを"使ってる"のがここ。KD-Treeを使うと、この部分が効率的に計算できる
    }

    // KD-Treeに関するオペレータ

    KDTreeNode* create_kdtree_sub(){
        // KD-Treeを"作ってる"のがここ
    }
public:
    // KD-Treeのコンストラクタ/デストラクタ

    // KD-Tree使用の発火場所
    void SearchKNN(typename KDTree::ResultQueue* pqueue, typename KDTree<T>::Query& query){
        locate_points(pqueue, root, query);
    }
    void AddPoint(const T& point){
        points.push_back(point);
    }

    // KD-Tree作成の発火場所
    void CreateKDTree(){
        root = create_kdtree_sub(points.begin(), points.end(), 0);
    }
};

inline bool intersect_scene(){
    // rayを渡してオブジェクトに当たるかどうかを返す
    // もし当たるならどのオブジェクトかを返す関数
}

typedef KDTree<Photon> PhotonMap;
void create_photon_map(const int shoot_photon_num, PhotonMap *photon_map){
    // フォトンマップの作成。ほぼレイトレ
    // この段階ではvectorにプッシュしてるだけ

    // 最後でvectorに詰めた情報をKD-Treeに整形
    photon_map->CreateKDTree();
    std::cout << "Done." << std::endl;
}

Color radiance(){
    // KD-Treeを使って輝度を計算する。SearchKNNとかが出てくるのがここ
}

// 画像を出力するための色々

int main(){

}