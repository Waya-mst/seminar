#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <queue>

const double PI = 3.14159265358979323846;
const double INF = 1e20;
const double EPS = 1e-6;
const double MaxDepth = 5;

// 返すのは0.0~1.0までの値
inline double rand01() {return (double)rand()/RAND_MAX;}

struct Vec {
    double x{};
    double y{};
    double z{};

    constexpr Vec() noexcept = default;
    constexpr Vec(
        double x_,
        double y_,
        double z_
    ) noexcept : x{x_}, y{y_}, z{z_} {}

    constexpr Vec operator+(const Vec& b) const noexcept{
        return {x + b.x, y + b.y, z + b.z};
    }

    constexpr Vec operator-(const Vec& b) const noexcept{
        return {x - b.x, y - b.y, z - b.z};
    }

    constexpr Vec operator*(const double& b) const noexcept{
        return {x * b, y * b, z * b};
    }

    constexpr Vec operator/(const double& b) const noexcept{
        return {x / b, y / b, z / b};
    }

    constexpr double LengthSquared() const {return x*x + y*y + z*z;}
    constexpr double Length() const {return sqrt(LengthSquared());}
};

inline Vec operator*(double f, const Vec &v){return f*v;}
inline Vec Normalize(const Vec &v){return v/v.Length();}

inline const Vec Multiply(const Vec &v1, const Vec &v2){
    return Vec(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z); 
}

inline const double Dot(const Vec &v1, const Vec &v2){
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline const Vec Cross(const Vec &v1, const Vec &v2){
    return Vec((v1.y * v2.z) - (v1.z * v2.y), (v1.z * v2.x) - (v1.x * v2.z), (v1.x * v2.y) - (v1.y * v2.x));
}

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

    Sphere(
        const double radius_,
        const Vec &pos_,
        const Color &emission_, const Color &color_,
        const ReflectionType ref_type_
    ): radius(radius_), pos(pos_), emission(emission_), color(color_), ref_type(ref_type_) {}

    const double intersect(const Ray &ray){
        Vec o_p = pos.operator-(ray.org);
        const double b = Dot(o_p, ray.dir);
        const double det = b*b - Dot(o_p, o_p) + radius * radius;
        if(det >= 0.0){
            const double sqrt_det = sqrt(det);
            const double t1 = b - sqrt_det;
            const double t2 = b + sqrt_det;
            if(t1 > EPS){
                return t1;
            }else{
                return t2;
            }
        }
        return 0.0;
    }
};

// シーンデータ
static std::vector<Sphere*> spheres = {
    new Sphere(5.0, Vec(50.0, 75.0, 81.6),Color(12,12,12), Color(), DIFFUSE),//照明
    new Sphere(1e5, Vec( 1e5+1,40.8,81.6), Color(), Color(0.75, 0.25, 0.25),DIFFUSE),// 左
    new Sphere(1e5, Vec(-1e5+99,40.8,81.6),Color(), Color(0.25, 0.25, 0.75),DIFFUSE),// 右
    new Sphere(1e5, Vec(50,40.8, 1e5),     Color(), Color(0.75, 0.75, 0.75),DIFFUSE),// 奥
    new Sphere(1e5, Vec(50,40.8,-1e5+170), Color(), Color(), DIFFUSE),// 手前
    new Sphere(1e5, Vec(50, 1e5, 81.6),    Color(), Color(0.75, 0.75, 0.75),DIFFUSE),// 床
    new Sphere(1e5, Vec(50,-1e5+81.6,81.6),Color(), Color(0.75, 0.75, 0.75),DIFFUSE),// 天井
    new Sphere(16.5,Vec(27,16.5,47),       Color(), Color(0,0,0.3), SPECULAR),// 鏡
    new Sphere(16.5,Vec(73,16.5,78),       Color(), Color(0,0.3,0.3), REFRACTION),//ガラス
};
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
    struct Query{
        Vec search_position; // 探索中心
        Vec normal; // 探索中心における法線
        double max_distance2; // 探索の最大半径
        size_t max_search_num; // 最大探索点数
        Query(
            const Vec& search_position_,
            const Vec& normal_,
            const double max_distance2_,
            const size_t max_search_num_
        ): search_position(search_position_),
        normal(normal_),
        max_distance2(max_distance2_),
        max_search_num(max_search_num_){}
    };

    struct ElementForQueue{
        const T* point;
        double distance2;
        
        ElementForQueue(const T* point_, double distance2_)
        : point(point_), distance2(distance2_){}

        bool operator<(const ElementForQueue& b) const {
            return distance2 < b.distance2;
        }
    };

    typedef std::priority_queue<ElementForQueue, std::vector<ElementForQueue>> ResultQueue;

private:
    std::vector<T> points;
    struct KDTreeNode{
        T* point; // フォトンの格納場所
        KDTreeNode* left; // leftへのポインタ
        KDTreeNode* right; // rightへのポインタ
        int axis; // x,y,z軸どれを基準に左右に分けたか
    };
    KDTreeNode* root;
    void delete_kdtree(KDTreeNode* node){
        if(node == NULL) return;
        delete_kdtree(node->left);
        delete_kdtree(node->right);
        delete node;
    }

    // KD-Treeの中からクエリの条件に合うフォトンを探索してResultQueueに詰める
    // 実質KD-Treeを"使ってる"のがここ。KD-Treeを使うと、この部分が効率的に計算できる
    void locate_points(typename KDTree<T>::ResultQueue* pqueue, KDTreeNode* node, typename KDTree<T>::Query &query){
        if(node == NULL) return;
        const int axis = node->axis;

        double delta;
        switch(axis){
            case 0: delta = query.search_position.x - node->point->pos.x; break;
            case 1: delta = query.search_position.y - node->point->pos.y; break;
            case 2: delta = query.search_position.z - node->point->pos.z; break;
        }

        const Vec dir = node->point->pos - query.search_position;
        const double distance2 = dir.LengthSquared();
        const double dt = Dot(query.normal, dir/sqrt(distance2));
        if(distance2 < query.max_distance2 && fabs(dt) <= query.max_distance2 * 0.01){
            // ここで探索結果をResultQueueに詰める
            pqueue->push(ElementForQueue(node->point, distance2));
            if(pqueue->size() > query.max_search_num){
                // 詰めすぎたら破棄
                pqueue->pop();
                // 探索範囲を狭める
                query.max_distance2 = pqueue->top().distance2;
            }
        }
        if(delta > 0.0){
            // 探索中心がノードの?側にあるので、ノードの?側を優先して探索
            locate_points(pqueue, node->/*right or left*/, query);

            // 分割平面と探索中心の距離が探索範囲より遠ければ、探索中心から見て分割平面の向こう側の点は調べても意味ない(絶対探索範囲より遠くなる)
            if(/*???*/){
                locate_points(pqueue, node->/*right or left*/, query);
            }
        }else{
            locate_points(pqueue, node->/*right or left*/, query);
            if(/*???*/){
                locate_points(pqueue, node->/*right or left*/, query);
            }
        }
    }

    static bool kdtree_less_operator_x(const T& left, const T& right){
        return left.pos.x < right.pos.x;
    }
    static bool kdtree_less_operator_y(const T& left, const T& right){
        return left.pos.y < right.pos.y;
    }
    static bool kdtree_less_operator_z(const T& left, const T& right){
        return left.pos.z < right.pos.z;
    }
    
    // KD-Treeを"作ってる"のがここ。
    KDTreeNode* create_kdtree_sub(typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end, int depth){
        if(end - begin <= 0){
            return NULL;
        }
        const int axis = depth % 3;

        switch(axis){
            case 0: std::sort(begin, end, kdtree_less_operator_x); break;
            case 1: std::sort(begin, end, kdtree_less_operator_y); break;
            case 2: std::sort(begin, end, kdtree_less_operator_z); break;
        }
        const int median = (end-begin) / 2;
        KDTreeNode* node = new KDTreeNode;
        node->axis = axis;
        node->point = &(*(begin+median));

        node->left = create_kdtree_sub(begin, begin+median, depth+1);
        node->right = create_kdtree_sub(begin+median+1, end, depth+1);
        return node;
    }
public:
    KDTree(){
        root = NULL;
    }
    virtual ~KDTree(){
        delete_kdtree(root);
    }
    size_t Size(){
        return points.size();
    }

    void SearchKNN(typename KDTree::ResultQueue* pqueue, typename KDTree<T>::Query& query){
        locate_points(pqueue, root, query);
    }
    void AddPoint(const T& point){
        points.push_back(point);
    }
    void CreateKDTree(){
        root = create_kdtree_sub(points.begin(), points.end(), 0);
    }
};

inline bool intersect_scene(const Ray& ray, double *t, int *id){
    const size_t n = spheres.size();
    *t = INF;
    *id = -1;
    for(int i = 0; i < int(n); i++){
        double d = spheres[i]->intersect(ray);
        if(d > 0.0 && d < *t){
            *t = d;
            *id = i;
        }
    }
    return *t < INF;
}

typedef KDTree<Photon> PhotonMap;
void create_photon_map(const int shoot_photon_num, PhotonMap *photon_map){
    std::cout << "Shooooooooting photons... (" << shoot_photon_num << " photons)" << std::endl;
    for(int i = 0; i < shoot_photon_num; i++){
        const double r1 = 2 * PI * rand01();
        const double r2 = 1.0 - 2.0 * rand01();
        
        const Vec light_pos = spheres[LightID]->pos +
        (Vec(sqrt(1.0 - r2*r2) * cos(r1), sqrt(1.0 - r2*r2) * sin(r1), r2) * (spheres[LightID]->radius + EPS));
        const Vec normal = Normalize(light_pos - spheres[LightID]->pos);
        Vec w, u, v;
        w = normal;
        if(fabs(w.x) > 0.1){
            u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
        }else{
            u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
        }
        v = Cross(w, u);
        
        const double u1 = 2*PI*rand01();
        const double u2 = rand01(), u2s = sqrt(u2);
        Vec light_dir = Normalize((u * cos(u1) * u2s + v * sin(u1) * u2s + w * sqrt(1.0 - u2)));

        Ray now_ray(light_pos, light_dir);
        Color now_flux = spheres[LightID]->emission * 4.0 * PI * pow(spheres[LightID]->radius, 2.0) / shoot_photon_num;
        
        bool trace_end = false;
        for(; !trace_end;){
                if(std::max(now_flux.x, std::max(now_flux.y, now_flux.z)) <= 0.0)
                    break;
            
            double t;
            int id;
            std::cout << "creating map" << std::endl;
            if(! intersect_scene(now_ray, &t, &id))
                break;
            const Sphere &obj = *spheres[id];
            const Vec hitpoint = now_ray.org + now_ray.dir * t;
            const Vec normal = Normalize(hitpoint - obj.pos);
            const Vec orienting_normal = Dot(normal, now_ray.dir) < 0.0 ? normal : (normal * -1.0);
            switch(obj.ref_type){
                case DIFFUSE:{
                    photon_map->AddPoint(Photon(hitpoint, now_flux, now_ray.dir));

                    const double probability = (obj.color.x + obj.color.y + obj.color.z) / 3;
                    if(probability > rand01()){
                        Vec w, u, v;
                        w = orienting_normal;
                        if(fabs(w.x) > 0.1)
                            u = Normalize(Cross(Vec(0.0, 1.0, 0.0), w));
                        else
                            u = Normalize(Cross(Vec(1.0, 0.0, 0.0), w));
                        v = Cross(w, u);
                        
                        const double r1 = 2 * PI * rand01();
                        const double r2 = rand01(), r2s = sqrt(r2);
                        Vec dir = Normalize((u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.0 - r2)));

                        now_ray = Ray(hitpoint, dir);
                        // 
                        now_flux = Multiply(now_flux, obj.color) / probability;
                        continue;
                    }else{
                        trace_end = true;
                        continue;
                    }
                }break;
                case SPECULAR:{
                    now_ray = /*反射するレイのベクトル*/;
                    now_flux = Multiply(now_flux, obj.color);
                    continue;
                }break;
                case REFRACTION: {
                    Ray reflection_ray = /*反射するレイのベクトル*/;
                    bool into = Dot(normal, orienting_normal) > 0.0;

                    const double nc = 1.0;
                    const double nt = 1.5;
                    const double nnt = into ? nc / nt : nt / nc;
                    const double ddn = Dot(now_ray.dir, orienting_normal);
                    const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                    if(cos2t < 0.0){// 全反射
                        now_ray = reflection_ray;
                        now_flux = Multiply(now_flux, obj.color);
                        continue;
                    }
                    Vec tdir = Normalize(
                        now_ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t))
                    );

                    const double a = nt - nc, b = nt + nc;
                    const double R0 = (a*a) / (b*b);
                    const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
                    const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
                    const double Tr = 1.0 - Re;
                    const double probability = Re;

                    if(rand01() < probability){
                        now_ray = reflection_ray;
                        now_flux = Multiply(now_flux, obj.color) * Re / probability;
                        continue;
                    }else{
                        now_ray = Ray(hitpoint, tdir);
                        now_flux = Multiply(now_flux, obj.color);
                        continue;
                    }
                }break;
            }
        }
    }
    std::cout << "Done. (" << photon_map->Size() <<  " photons are stored)" << std::endl;
    std::cout << "Creating KD-tree..." << std::endl;
    photon_map->CreateKDTree();
    std::cout << "Done." << std::endl;
}

Color radiance(const Ray& ray,
    const int depth,
    PhotonMap* photon_map,
    const double gather_radius,
    const int gather_max_photon_num)
{
        double t;
        int id;
        if(!intersect_scene(ray, &t, &id))
            return BackGroundColor;
        const Sphere &obj = *spheres[id];
        const Vec hitpoint = ray.org + ray.dir * t;
        const Vec normal = Normalize(hitpoint - obj.pos);
        const Vec orienting_normal = Dot(normal, ray.dir) > 0 ? normal : (normal * -1.0);
        double russian_roulette_probability = std::max(obj.color.x, std::max(obj.color.y, obj.color.z));
        if(depth > MaxDepth){
            if(rand01() >= russian_roulette_probability)
                return obj.emission;
        }else
            russian_roulette_probability = 1.0;

        switch(obj.ref_type){
            case DIFFUSE:{
                PhotonMap::ResultQueue pqueue;
                PhotonMap::Query query(hitpoint, orienting_normal, gather_radius, gather_max_photon_num);
                photon_map->SearchKNN(&pqueue, query);
                Color accumlated_flux;
                double max_distance2 = -1;

                std::vector<PhotonMap::ElementForQueue> photons;
                photons.reserve(pqueue.size());
                for(; !pqueue.empty();){
                    PhotonMap::ElementForQueue p = pqueue.top(); pqueue.pop();
                    photons.push_back(p);
                    // gather_max_photon_num個のうち,最も離れた点までの距離を使って輝度を計算
                    max_distance2 = std::max(max_distance2, p.distance2);
                }
                const double max_distance = sqrt(max_distance2);
                const double k = 1.1;
                for(int i = 0; i < photons.size(); i++){
                    //　ここの円錐フィルタってなに？
                    const double w = 1.0 - (sqrt(photons[i].distance2) / (k * max_distance));
                    const Color v = Multiply(obj.color, photons[i].point->power) / PI;
                    accumlated_flux = accumlated_flux + v * w;
                }
                accumlated_flux = accumlated_flux / (1.0 - 2.0 / (3.0 * k));
                if(max_distance2 > 0.0){
                    return obj.emission + accumlated_flux / (PI * max_distance2) / russian_roulette_probability;
                    // なんでrussian_roulette_probabilityで割ってるの？
                }
            }break;
            
            case SPECULAR: {
                // return obj.emission + 
                // radiance(Ray(hitpoint, ray.dir + orienting_normal * Dot(ray.dir, orienting_normal) * 2.0),
                //         depth + 1, photon_map, gather_radius, gather_max_photon_num) / russian_roulette_probability;
                return obj.emission +
                radiance(Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir)), depth+1, photon_map, gather_radius, gather_max_photon_num);
            }break;

            case REFRACTION: {
                Ray reflection_ray = Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir));
                bool into = Dot(normal, orienting_normal) > 0.0;

                const double nc = 1.0; // 真空の屈折率
                const double nt = 1.5; // オブジェクトの屈折率
                const double nnt = into ? nc / nt : nt / nc;
                const double ddn = Dot(ray.dir, orienting_normal);
                const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

                if(cos2t < 0.0){
                    return obj.emission + 
                    Multiply(obj.color, radiance(Ray(hitpoint, ray.dir - normal * 2.0 * Dot(normal, ray.dir)), depth + 1, photon_map, gather_radius, gather_max_photon_num));
                }
                Vec tdir = Normalize(ray.dir * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));

                // SchlickによるFresnelの反射係数の近似
                const double a = nt - nc, b = nt + nc;
                const double R0 = (a * a) / (b * b);
                const double c = 1.0 - (into ? -ddn : Dot(tdir, normal));
                const double Re = R0 + (1.0 - R0) * pow(c, 5.0);
                const double Tr = 1.0 - Re; // 屈折光の運ぶ光の量
                const double probability  = 0.25 + 0.5 * Re;

                if(depth > 2){
                    if(rand01() < probability){// 反射
                        return obj.emission +
                        Multiply(obj.color, radiance(reflection_ray, depth + 1, photon_map, gather_radius, gather_max_photon_num) * Re)
                        / probability
                        / russian_roulette_probability;
                    }else{// 屈折
                        return obj.emission +
                        Multiply(obj.color, radiance(reflection_ray, depth+1, photon_map, gather_radius, gather_max_photon_num) * Tr)
                        / (1.0 - probability)
                        / russian_roulette_probability;
                    }
                }else{
                    return obj.emission +
                    Multiply(obj.color, radiance(reflection_ray, depth+1, photon_map, gather_radius, gather_max_photon_num) * Re
                    + radiance(Ray(hitpoint, tdir), depth+1, photon_map, gather_radius, gather_max_photon_num) * Tr / russian_roulette_probability);
                }
            }break;
        }
        return Color();
}

struct HDRPixel {
    unsigned char r, g, b, e;
    HDRPixel(
        const unsigned char r_ = 0,
        const unsigned char g_ = 0,
        const unsigned char b_ = 0,
        const unsigned char e_ = 0
    ) :
    r(r_), g(g_), b(b_), e(e_) {};

    unsigned char get(int idx) {
        switch (idx) {
            case 0: return r;
            case 1: return g;
            case 2: return b;
            case 3: return e;
        } return 0;
    }
};

    // doubleのRGB要素を.hdrフォーマット用に変換
HDRPixel get_hdr_pixel(const Color &color) {
  double d = std::max(color.x, std::max(color.y, color.z));
  if (d <= 1e-32)
    return HDRPixel();
  int e;
  double m = frexp(d, &e); // d = m * 2^e
  d = m * 256.0 / d;
  return HDRPixel(color.x * d, color.y * d, color.z * d, e + 128);
}

// 書き出し用関数
void save_hdr_file(const std::string &filename, const Color* image, const int width, const int height) {
  FILE *fp = fopen(filename.c_str(), "wb");
  if (fp == NULL) {
    std::cerr << "Error: " << filename << std::endl;
    return;
  }
  // .hdrフォーマットに従ってデータを書きだす
  // ヘッダ
  unsigned char ret = 0x0a;
  fprintf(fp, "#?RADIANCE%c", (unsigned char)ret);
  fprintf(fp, "# Made with 100%% pure HDR Shop%c", ret);
  fprintf(fp, "FORMAT=32-bit_rle_rgbe%c", ret);
  fprintf(fp, "EXPOSURE=1.0000000000000%c%c", ret, ret);

  // 輝度値書き出し
  fprintf(fp, "-Y %d +X %d%c", height, width, ret);
  for (int i = height - 1; i >= 0; i --) {
    std::vector<HDRPixel> line;
    for (int j = 0; j < width; j ++) {
      HDRPixel p = get_hdr_pixel(image[j + i * width]);
      line.push_back(p);
    }
    fprintf(fp, "%c%c", 0x02, 0x02);
    fprintf(fp, "%c%c", (width >> 8) & 0xFF, width & 0xFF);
    for (int i = 0; i < 4; i ++) {
      for (int cursor = 0; cursor < width;) {
	const int cursor_move = std::min(127, width - cursor);
	fprintf(fp, "%c", cursor_move);
	for (int j = cursor;  j < cursor + cursor_move; j ++)
	  fprintf(fp, "%c", line[j].get(i));
	cursor += cursor_move;
      }
    }
  }

  fclose(fp);
}

int main(int argc, char **argv){
    int width = 640;
    int height = 480;
    int photon_num = 50000;
    double gather_photon_radius = 32.0;
    int gather_max_photon_num = 3;

    Ray camera(Vec(50.0, 52.0, 295.6), Normalize(Vec(0.0, -0.04, -1.0)));
    Vec cx = Vec(width * 0.5 / height, 0.0, 0.0);
    Vec cy = Normalize(Cross(cx, camera.dir)) * 0.5;
    Color *image = new Color[width * height];

    PhotonMap photon_map;
    create_photon_map(photon_num, &photon_map);

    for(int y = 0; y < height; y++){
        std::cerr << "Rendering " << (100.0 * y / (height - 1)) << "%" << std::endl;
        srand(y * y * y);
        
        for(int x = 0; x < width; x++){
            int image_index = y * width + x;
            image[image_index] = Color();
            for(int sy = 0; sy < 2; sy++){
                for(int sx = 0; sx < 2; sx++){
                    const double r1 = 2.0 * rand01(), dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
                    const double r2 = 2.0 * rand01(), dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);
                    Vec dir = cx * (((sx + 0.5 + dx) / 2.0 + x) / width - 0.5) +
                    cy * (((sy + 0.5 + dy) / 2.0 + y) / height - 0.5) + camera.dir;
                    Color accum = radiance(Ray(camera.org + dir * 130.0, Normalize(dir)), 3, &photon_map, gather_photon_radius, gather_max_photon_num);
                    image[image_index] = image[image_index] + accum;
                }
            }
        }
    }
    save_hdr_file(std::string("image.hdr"), image, width, height);
    return 0;
}


