#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "gen_points.h"

typedef struct tree_node {
    int center_idx;
    double radius;
    long left;
    long right;
} node;

double **centers;

int n_dims; // Dimensions of the input points
long np;    // Number of generated points
/* Array of points: allocated at the beggining of the program;
* Points are randomly generated at the start of the program and don't change as it runs
*/
double **pt_array;
/* Struct that saves all tree information, for later printing;
* Has 2 * np - 1 points
*/
node *tree;
double *proj_scalar;
/* vector with the sorted indices. 
* All the computations are made in this vector:
* sorting according to projection value, splitting in right/left ball
* ...
*/
long *idx;

/* Function to print a point to a given file */
void print_point(double *pt, int dim, FILE *fd) {
    for (int i = 0; i < dim; ++i)
        fprintf(fd, "%f ", pt[i]);
    fprintf(fd, "\n");
}

/* Computes the squared distance between 2 points */
double dist(double *pt1, double *pt2) {
    double dist = 0.0;
    for (int d = 0; d < n_dims; ++d)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return dist;
}

void swap(long *i, long *j) {
    long temp = *i;
    *i = *j;
    *j = temp;
}

void furthest_apart(int l, int r, long *a_idx, long *b_idx) {
    double max_d = 0.0;
    double d;
    long idx0 = idx[l];
    for (int i = l + 1; i < r; ++i) {
        if (idx[i] < idx0)
            idx0 = idx[i];
    }
    for (int i = l; i < r; ++i) {
        d = dist(pt_array[idx0], pt_array[idx[i]]);
        if (d > max_d) {
            max_d = d;
            *a_idx = idx[i];
        }
    }
    max_d = 0;
    for (int i = l; i < r; ++i) {
        d = dist(pt_array[*a_idx], pt_array[idx[i]]);
        if (d > max_d) {
            max_d = d;
            *b_idx = idx[i];
        }
    }
    if (pt_array[*a_idx][0] > pt_array[*b_idx][0])
        swap(a_idx, b_idx);
}

double orth_projv1(double *a, double *b, double *p) {
    // Orthogonal projection in the first coordinate
    double proj = 0.0;
    for (int i = 0; i < n_dims; ++i)
        proj += (p[i] - a[i]) * (b[i] - a[i]);

    return proj;
}

void orth_projv2(double *a, double *b, long idx_p, double *ab_proj) {
    double abnorm = 0.0, aux, u = proj_scalar[idx_p];
    for (int i = 0; i < n_dims; ++i) {
        aux = (b[i] - a[i]);
        ab_proj[i] = u * aux;
        abnorm += aux * aux;
    }
    for (int i = 0; i < n_dims; ++i)
        ab_proj[i] = ab_proj[i] / abnorm + a[i];
}

int is_seq(int i, int j) {
    return proj_scalar[i] <= proj_scalar[j];
}

/* Gets the k-th element of a piece of the projections array (from l to h):
* changes only the idx (indices) array;
* the comparison critirion is the first coordinate of all projections (stored in proj_scalar);
* also partitions the data: at the left of the k-th element of the indices array are all the indices
*   corresponding to points whose the first coordinate projection in the current ab line is not bigger
*   than that of the k-th element (points that should be on the left ball);
* at the right of the k-th element of the indices array are all the indices
*   corresponding to points whose the first coordinate projection in the current ab line is not smaller
*   than that of the k-th element (points that should be on the right ball).
*/
int get_kth_element(int k, int l, int h) {

    if (h == l + 1)
        return idx[l];

    else {
        int i = l;
        int j = h;

        while (i < j) {
            do {
                i++;
            } while (i < h && is_seq(idx[i], idx[l]));
            if (i == h) {
                j = h - 1;
                break;
            }
            do {
                j--;
            } while (!is_seq(idx[j], idx[l]));
            if (i < j)
                swap(&idx[i], &idx[j]);
        }

        swap(&idx[l], &idx[j]);

        if (j - l == k)
            return idx[j];
        if (k < j - l)
            return get_kth_element(k, l, j);
        else
            return get_kth_element(k - (j - l) - 1, j + 1, h);
    }
}

void get_median(int l, int h, int *median_1, int *median_2) {
    if ((h - l) % 2 == 1)
        *median_1 = get_kth_element((h - l - 1) / 2, l, h);
    else {
        *median_1 = get_kth_element((h - l) / 2 - 1, l, h);
        *median_2 = get_kth_element(0, l + (h - l) / 2, h);
    }
}

double *abproj;

long n_nodes = 0; // Total number of tree nodes (in the end will equal 2 * np - 1)
long n_center = 0;

long ballAlg(long l, long r) {

    ++n_nodes;
    long id = n_nodes - 1;

    if (r - l == 1) {
        tree[id].center_idx = idx[l];
        tree[id].radius = 0;
        tree[id].left = -1;
        tree[id].right = -1;
        return id;
    }

    long c_id = n_center++;

    long a, b;
    // 2. Compute points a and b, furthest apart in the current set (approx)
    furthest_apart(l, r, &a, &b);

    // 3. Perform the orthogonal projection of all points onto line ab
    for (int i = l; i < r; ++i)
        proj_scalar[idx[i]] = orth_projv1(pt_array[a], pt_array[b], pt_array[idx[i]]);

    // 4. Compute the center, defined as the median point over all projections
    int m1, m2 = -1;
    get_median(l, r, &m1, &m2);

    if ((r - l) % 2) {
        orth_projv2(pt_array[a], pt_array[b], m1, centers[c_id]);
    } else {
        double abnorm = 0.0, aux, u = (proj_scalar[m1] + proj_scalar[m2]) / 2;
        for (int i = 0; i < n_dims; ++i) {
            aux = (pt_array[b][i] - pt_array[a][i]);
            centers[c_id][i] = u * aux;
            abnorm += aux * aux;
        }
        for (int i = 0; i < n_dims; ++i)
            centers[c_id][i] = centers[c_id][i] / abnorm + pt_array[a][i];
    }

    double max_r = 0, rad;
    for (int i = l; i < r; ++i) {
        rad = dist(centers[c_id], pt_array[idx[i]]);
        if (rad > max_r)
            max_r = rad;
    }
    tree[id].center_idx = c_id;
    tree[id].radius = sqrt(max_r);

    tree[id].left = ballAlg(l, l + (r - l) / 2);
    tree[id].right = ballAlg(l + (r - l) / 2, r);

    return id;
}

void print_tree(node *tree) {
    fprintf(stdout, "%d %ld\n", n_dims, 2 * np - 1);
    for (long i = 0; i < 2 * np - 1; ++i) {
        fprintf(stdout, "%ld %ld %ld %f ", i, tree[i].left, tree[i].right, tree[i].radius);
        if (tree[i].left == -1)
            print_point(pt_array[tree[i].center_idx], n_dims, stdout);
        else
            print_point(centers[tree[i].center_idx], n_dims, stdout);
    }
}

int main(int argc, char **argv) {

    double exec_time;
    exec_time = -omp_get_wtime();

    // 1. Get input sample points
    pt_array = get_points(argc, argv, &n_dims, &np);

    idx = (long *)malloc(sizeof(long) * np);
    for (int i = 0; i < np; ++i)
        idx[i] = i;

    proj_scalar = (double *)malloc(np * sizeof(double));

    tree = (node *)malloc((2 * np - 1) * sizeof(node));
    double *center_aux = (double *)malloc((np - 1) * n_dims * sizeof(double));
    centers = (double **)malloc((np - 1) * sizeof(double *));
    for (int i = 0; i < np - 1; ++i)
        centers[i] = &center_aux[i * n_dims];

    ballAlg(0, np);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1lf\n", exec_time);

    print_tree(tree);

    free(pt_array[0]);
    free(pt_array);

    free(centers[0]);
    free(centers);

    free(tree);

    free(proj_scalar);
    free(idx);

    exit(0);
}
