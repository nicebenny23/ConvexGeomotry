// ConvexGeomotry.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <vector>
#include <iostream>
#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <optional>
#include "random.h"
#include <unordered_set>
struct Vector {
    std::vector<double> coords;

    Vector() = default;

    explicit Vector(const std::vector<double>& v)
        : coords(v) {
    }
    explicit Vector(std::vector<double>&& v)
        : coords(std::move(v)) {
    }

    explicit Vector(size_t n)
        : coords(n, 0.0) {
    }

    size_t dim() const {
        return coords.size();
    }

    static void check_dimensions(const Vector& a, const Vector& b) {
        if (a.dim() != b.dim())
            throw std::invalid_argument("Point dimensions do not match.");
    }

    static double dot_product(const Vector& a, const Vector& b) {
        check_dimensions(a, b);

        double sum = 0.0;
        for (size_t i = 0; i < a.dim(); i++)
            sum += a[i] * b[i];

        return sum;
    }

    double magnitude() const {
        return std::sqrt(dot_product(*this, *this));
    }
    Vector normal() const {
        return *this / magnitude();
    }
    Vector& normalize() {
        *this /= magnitude();
        return *this;
    }
    double& operator[](size_t i) {
        return coords[i];
    }

    const double& operator[](size_t i) const {
        return coords[i];
    }
    
    Vector(const Vector&) = default;

    Vector(Vector&&) = default;
    Vector operator+(const Vector& other) const {
        check_dimensions(*this, other);

        Vector result(dim());
        for (size_t i = 0; i < dim(); i++)
            result[i] = coords[i] + other[i];

        return result;
    }
 
    Vector& operator=(const Vector& other) {
        check_dimensions(*this, other);
        coords = other.coords;
        return *this;
    }
    void abs_each() {
        for (size_t i = 0; i < dim(); i++)
        {
            coords[i] = std::abs(coords[i]);
        }
    }
    Vector& operator=(Vector&& other) {
        check_dimensions(*this, other);
        coords = std::move(other.coords);
        return *this;
    }
    Vector& operator+=(const Vector& other) {
        return *this = *this + other;
    }

    Vector operator-() const {
        return *this * -1.0;
    }

    Vector operator-(const Vector& other) const {
        return *this + (-other);
    }

    Vector& operator-=(const Vector& other) {
        return *this = *this - other;
    }

    Vector operator*(double scalar) const {
        Vector result(dim());

        for (size_t i = 0; i < dim(); i++)
            result[i] = coords[i] * scalar;

        return result;
    }

    Vector& operator*=(double scalar) {
        for (size_t i = 0; i < dim(); i++)
            coords[i] *= scalar;

        return *this;
    }

    Vector operator/(double scalar) const {
        if (scalar == 0.0)
            throw std::invalid_argument("Division by zero.");

        return *this * (1.0 / scalar);
    }

    Vector& operator/=(double scalar) {
        return *this *= 1/ scalar;
    }

    bool operator==(const Vector& other) const {
        return coords == other.coords;
    }

    bool operator!=(const Vector& other) const {
        return !(*this == other);
    }
    static Vector random(size_t dim) {
        std::vector<double> list;
        for (size_t i = 0; i < dim; i++)
        {
            list.push_back(gaussian());
        }
        Vector res = Vector(std::move(list));
        res.normalize();
        return res;
    }

    Vector project_onto_this(const Vector& projected_onto) const {
        return normal()*Vector::dot_product(projected_onto,normal());
    }
};namespace std {
    template<>
    struct hash<Vector> {
        size_t operator()(const Vector& v) const noexcept {
            size_t seed = 0;

            for (double x : v.coords) {
                size_t h = std::hash<double>{}(x);

                // hash combine (boost-style)
                seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }

            return seed;
        }
    };
#include <ostream>

        std::ostream & operator<<(std::ostream & os, const Vector & v)
    {
        os << '(';

        for (size_t i = 0; i < v.dim(); ++i)
        {
            if (i != 0)
                os << ", ";

            os << v[i];
        }

        os << ')';
        return os;
    }
}

struct HalfSpace {
    Vector normal;
    double bound;
    HalfSpace(const Vector& norm, double bnd) :normal(norm), bound(bnd) {
        bound /= norm.magnitude();
        normal.normalize();
    }
    HalfSpace operator-() {
        return HalfSpace(-normal, -bound);
    }
    bool contains(const Vector& other) const {
        return Vector::dot_product(normal, other) <= bound;
    }
    bool suffieciently_close(const Vector& other) const {
        return std::abs(Vector::dot_product(normal, other)- bound)<= 1e-6;
    }
};
inline std::vector<Vector> gramm_shmitt_process(std::vector<Vector> vectors) {
    for (int i =0;i<vectors.size();i++)
    {
        for (size_t j = 0; j<i; j++)
        {
            vectors[i]-=vectors[j].project_onto_this(vectors[i]);
        }
    }
    return vectors;
}
    inline HalfSpace space_from_vectors(std::vector<Vector> vectors) {
        std::vector<Vector> edges;

        for (size_t i = 1; i < vectors.size(); i++)
            edges.push_back(vectors[i] - vectors[0]);

        auto basis = gramm_shmitt_process(edges);

        Vector normal = Vector::random(vectors[0].dim());

        for (auto& b : basis)
        {
            normal -= b.project_onto_this(normal);
        }
        return HalfSpace(normal,Vector::dot_product(vectors[0], normal));
    }
    void space_from_vector_check(std::vector<Vector> vectors) {
        HalfSpace p1 = space_from_vectors(vectors);
        HalfSpace p2 = space_from_vectors(vectors);
        int l = 3;
    }
struct Frep {
    Frep() {

    }
    std::vector<HalfSpace> planes;
    void add(const HalfSpace& space) {
        planes.push_back(space);
    }
    bool contains(const Vector& point) const{
        for (const HalfSpace& h:planes)
        {
            if (!h.contains(point))
            {
                return false;
            }
        }
        return true;
    }
};
struct Simplex {
    std::vector<Vector> points;
    size_t dim() {
        return points[0].dim();
    }
    double volume() {
        std::vector<Vector> perp;
        for (size_t i = 1; i <= dim(); i++)
        {
            perp.push_back(points[i] - points[0]);
        }
        perp =gramm_shmitt_process(perp);
        double amt = 1;
        for (size_t i = 0; i <dim(); i++)
        {
            amt *= perp[i].magnitude()/(i+1);
        }
        return amt;
    }

};
inline void combinations_helper(
    int n,
    int k,
    int start,
    std::vector<size_t>& current,
    std::vector<std::vector<size_t>>& result)
{
    if (current.size() == k) {
        result.push_back(current);
        return;
    }

    for (int i = start; i <= n - (k - current.size()); ++i) {
        current.push_back(i);
        combinations_helper(n, k, i + 1, current, result);
        current.pop_back();
    }
}

inline std::vector<std::vector<size_t>> combinations(int n, int k)
{
    std::vector<std::vector<size_t>> result;
    std::vector<size_t> current;
    combinations_helper(n, k, 0, current, result);
    return result;
}
struct Vrep {
    size_t dim() const {
        return points[0].dim();
    }
    void add(Vector point) {
        points.push_back(point);
    }
    std::vector<Vector> points;
    std::pair<Vrep, Vrep> passing_failing(HalfSpace space) {
        Vrep passes;
        Vrep fails;
        HalfSpace neg_space = -space;
        for (const Vector& v : points)
        {
            if (space.suffieciently_close(v))
            {
                passes.add(v);
                fails.add(v);
            }
            else {
                if (space.contains(v))
                {
                    passes.add(v);
                }
                if (neg_space.contains(v))
                {
                    fails.add(v);
                }
            }
        }
        return std::pair<Vrep, Vrep>(passes, fails);
    }
    bool is_supporting_hyperplane(HalfSpace space) {
        size_t passes=0;
        size_t fails=0;
        HalfSpace neg_space = -space;
        for (const Vector& v : points)
        {
            if (!space.suffieciently_close(v))
            {
                if (space.contains(v))
                {
                    passes++;
                }
                else
                {
                    fails++;
                }
            }

            if (fails != 0 && passes != 0)
            {
                return false;
            }
        }
        if (size()-passes-fails!=dim())
        {
            int l = 3;
        }
        return true;
    }
    void hull() {
            std::unordered_set<Vector > point_list;
            std::vector<std::vector<size_t>> combines{ combinations(points.size(), dim() ) };
            for (const std::vector<size_t>& indices: combines)
            {
                std::vector<Vector> vertices;
                for (size_t i : indices) {
                    vertices.push_back(points[i]);
                }
                HalfSpace space = space_from_vectors(vertices);
                if (is_supporting_hyperplane(space))
                {
                    for (size_t i = 0; i <vertices.size(); i++)
                    {
                        if (!space.suffieciently_close(vertices[i]))
                        {
                            int l = 3;
                        }
                        point_list.emplace(vertices[i]);
                    }
                }
            }

            points.clear();
            for (const auto& x : point_list) {
                points.push_back(x);
            }
            
            int l = 3;
    }
    
    size_t size() const {
        return points.size();
    }
    Vector point_in() {
        Vector total(dim());
        for (size_t i = 0; i < size(); i++)
        {
            total += points[i];
        }
        return total / size();
    }
    std::vector<Simplex> simplexes() {
        std::vector<Simplex> point_list;
        std::vector<std::vector<size_t>> combines{ combinations(points.size(), dim()) };
        Vector center = point_in();
        for (const std::vector<size_t>& indices : combines)
        {
            if (indices[0]==5&&indices[1]==13)
            {
                int l = 4;
            }
            std::vector<Vector> vertices;
            for (size_t i : indices) {
                vertices.push_back(points[i]);
            }
            HalfSpace space = space_from_vectors(vertices);
            if (is_supporting_hyperplane(space))
            {

                    vertices.push_back(center); 
                    point_list.push_back(Simplex{ vertices });
            }
        }
        return point_list;
    }
    double volume() {
        hull();
        std::vector<Simplex> simplex_list= simplexes();
        double volume = 0;
        for (size_t i = 0; i < simplex_list.size(); i++)
        {
            volume += simplex_list[i].volume();
        }
        return volume;
    }
};
struct Box {
    Vector scale;
    Vector center;
    size_t dim() const {
        return scale.dim();
    }
    Vector random_pnt() const {
        std::vector<double> points;
        for (size_t i = 0; i < dim(); i++)
        {
            points.push_back((random() - .5) * scale[i] + center[i]);
        }
        return Vector(points);
    }
    static Box from_min_max(Vector min, Vector max) {
        Vector sub = max - min;
        sub.abs_each();
        return Box{ .scale =sub,.center = (max + min) / 2 };
    }
    double volume() const {
        double val = 1;
        for (size_t i = 0; i < dim(); i++)
        {
            val *= scale[i];
        }
        return val;
    }
};
struct SupportFunction {
    virtual size_t dim() const = 0;
    virtual double operator()(Vector amt) const = 0;
    HalfSpace half_space_from(Vector amt) const {
        return HalfSpace(amt, this->operator()(amt));
    }
    Box box() const {
        Vector jank(dim());
        Vector min(dim());
        Vector max(dim());
        for (size_t i = 0; i < dim(); i++)
        {
            jank[i] = 1;
            if (i!=0)
            {
                jank[i - 1] = 0;
            }
            min[i] = -operator()(-jank);
            max[i] = operator()(jank);
        }
        return Box::from_min_max(min, max);
    }
    Frep f_rep(size_t size) const {
        Frep p;
        for (size_t i = 0; i < size; i++)
        {
            Vector next = Vector::random(dim());

            p.add(half_space_from(next));
        }
        return p;
    }
    bool contains(Vector pnt, size_t evals) const {
        for (size_t i = 0; i < evals; i++)
        {
            Vector next = Vector::random(dim());
            if (!half_space_from(next).contains(pnt))
            {
                return false;
            }
        }
        return true;
    }

    double volume(size_t evals,size_t cont) const {
        Box eval_box = box();
        size_t count = 0;
        Frep rep = f_rep(cont);
        Vrep points;
        for (size_t i = 0; i < evals; i++)
        {
            Vector point = eval_box.random_pnt();
            if (rep.contains(point))
            {
                points.add(point);

            }
        }
        return points.volume();
    }

};
struct Sphere:SupportFunction{
    size_t dimention = 2;
    size_t dim() const {
        return dimention;
    }

    double operator()(Vector amt) const {
        return 1;
    }
};
int main()
{
    std::cout << Sphere().volume(1000,1000);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
