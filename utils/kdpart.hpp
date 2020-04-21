/**
* @file kdpart.hpp
* @brief kdtree based partitioning classes.
* @author Aurel Neic
* @version
* @date 2017-02-14
*/

#ifndef _KDPART_HPP
#define _KDPART_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#ifdef KDPART_MPI
#include <mpi.h>
#endif

namespace kdpart {

/// minimalistic internal point struct
template<class S>
struct vec3 {
  S x = S(), y = S(), z = S();

  void get(const S* p) {
    x = p[0], y = p[1], z = p[2];
  }
  void set(S* p) {
    p[0] = x, p[1] = y, p[2] = z;
  }
};

/// Bounding box struct
template<class S>
struct bbox {
  /// the bounds. bounds[0] = lower left (min), bounds[1] = upper right (max)
  kdpart::vec3<S> bounds[2];
};

/// split axis
enum axis {
  X = 0, Y, Z, UNSET
};

/// element definition
template<class T, class S>
struct elem {
  vec3<S> ctr;        ///< element center location
  T       eidx = -1;  ///< element index
};

/// the struct holding all partition data
template<class T, class S>
struct partition {
  T pidx = T(-1);                               ///< partition index
  T cnt  = T(-1);                               ///< (global) partition size
  kdpart::bbox<S> box;                          ///< (global) partition bounding box
  std::vector<elem<T,S> > * elems = NULL;       ///< (loval) elements in partition

  ~partition() {
    if(elems) delete elems;
  }
};

/**
* @brief Check if a given value is an integer power of two
*
*/
inline bool is_power_of_two(double val)
{
  double log_two     = log(val) / log(2.0);
  int    log_two_int = log_two;

  return (fabs(double(log_two_int) - log_two) < 1e-6);
}



/**
* @brief Combined floating ponit and integer pair
*
* Used by sequential_partitioner
*
*/
template<class T, class S>
struct mixed_pair {
  T v1;
  S v2;
};

/// sorting operator
template<class T, class S>
bool operator< (const mixed_pair<T,S> & lhs, const mixed_pair<T,S> & rhs)
{
  return lhs.v1 < rhs.v1;
}

/**
* @brief Clamp a value into an interval [start, end]
*
* @tparam V   Value type
* @tparam W   Interval boundary type
* @param val    The value we clamp
* @param start  The interval start value
* @param end    The interval end value
*
* @return The clamped value
*/
template<typename V, typename W>
V clamp(const V val, const W start, const W end) {
  if(val < start) return start;
  if(val > end) return end;
  return val;
}

/// Compute displacements from counts.
template<class T>
inline void dsp_from_cnt(const std::vector<T> & cnt, std::vector<T> & dsp)
{
  dsp.resize(cnt.size()+1);
  dsp[0] = 0;
  for(size_t i=0; i<cnt.size(); i++) dsp[i+1] = dsp[i] + cnt[i];
}

/// Compute counts from displacements.
template<class T>
inline void cnt_from_dsp(const std::vector<T> & dsp, std::vector<T> & cnt)
{
  cnt.resize(dsp.size() - 1);
  for(size_t i=0; i<dsp.size()-1; i++) cnt[i] = dsp[i+1] - dsp[i];
}

template<class V, class W>
void sort_copy(std::vector<V> & v1, std::vector<W> & v2)
{
  assert(v1.size() == v2.size());

  std::vector<mixed_pair<V,W> > pair_array(v1.size());

  for(size_t i=0; i<v1.size(); i++) {
    pair_array[i].v1 = v1[i];
    pair_array[i].v2 = v2[i];
  }

  std::sort(pair_array.begin(), pair_array.end());

  for(size_t i=0; i<v1.size(); i++) {
    v1[i] = pair_array[i].v1;
    v2[i] = pair_array[i].v2;
  }
}


/**
* @brief The value we add to the power of two if we need higher partitioning resolution.
*
* In the case that the number of processes (thus the number of requested partitions) is
* not a power of two, we partition with a higher partition number and then evenly
* distribute the partitions onto the processes. This reduces the difference between
* partition sizes.
*
*/
#define KD_ORDER_INC 5.0

#define KD_MIN_SIZE 64

#ifdef KDPART_MPI
template<class T, class S>
class parallel_partitioner
{
  private:
  // the kdtree_partitioner members
  std::vector<kdpart::partition<T,S> > _layout;  ///< the partitioning layout
  T _curpart;                              ///< current number of assigned partitions
  MPI_Comm _comm;                          ///< The used communicator, we get it from the mesh

  /// compute bounding box for a vector of elements
  inline kdpart::bbox<S> get_bbox(const std::vector<elem<T,S> > & elems)
  {
    double minmax[6], minmax_red[6];

    minmax[0] = 1e100;
    minmax[1] = 1e100;
    minmax[2] = 1e100;
    minmax[3] = -1e100;
    minmax[4] = -1e100;
    minmax[5] = -1e100;

    for(size_t i=0; i<elems.size(); i++) {
      kdpart::vec3<S> p = elems[i].ctr;

      if(minmax[0] > p.x) minmax[0] = p.x;
      if(minmax[1] > p.y) minmax[1] = p.y;
      if(minmax[2] > p.z) minmax[2] = p.z;
      if(minmax[3] < p.x) minmax[3] = p.x;
      if(minmax[4] < p.y) minmax[4] = p.y;
      if(minmax[5] < p.z) minmax[5] = p.z;
    }

    MPI_Allreduce(minmax, minmax_red, 3, MPI_DOUBLE, MPI_MIN, _comm);
    MPI_Allreduce(minmax+3, minmax_red+3, 3, MPI_DOUBLE, MPI_MAX, _comm);

    kdpart::bbox<S> box;
    kdpart::vec3<S> & min = box.bounds[0];
    kdpart::vec3<S> & max = box.bounds[1];

    min.x = minmax_red[0];
    min.y = minmax_red[1];
    min.z = minmax_red[2];
    max.x = minmax_red[3];
    max.y = minmax_red[4];
    max.z = minmax_red[5];

    return box;
  }

  /// get the longest axis of a bounding box
  inline kdpart::axis get_longest_axis(const kdpart::bbox<S> & box)
  {
    const kdpart::vec3<S> & min = box.bounds[0];
    const kdpart::vec3<S> & max = box.bounds[1];

    S x = fabs(min.x - max.x);
    S y = fabs(min.y - max.y);
    S z = fabs(min.z - max.z);

    return x > y && x > z ? X : y > z ? Y : Z;
  }

  /**
  * @brief Print some details on the layout after a splitting.
  *
  * The location of the last split will be highlighted.
  *
  * @param split_pos The partition where the last splitting occured.
  */
  inline void print_layout(const T split_pos = -1)
  {
    int rank; MPI_Comm_rank(_comm, &rank);

    if(rank != 0) return;

    if(split_pos > -1) {
      T idx = 0;
      while(idx < split_pos) {
        printf("%d : %d \n", int(_layout[idx].pidx), int(_layout[idx].cnt));
        idx++;
      }

      printf("----\n");
      printf("%d : %d \n", int(_layout[idx].pidx), int(_layout[idx].cnt));
      printf("%d : %d \n", int(_layout[idx+1].pidx), int(_layout[idx+1].cnt));
      printf("----\n");
      idx += 2;

      while(size_t(idx) < _layout.size() && _layout[idx].pidx > -1) {
        printf("%d : %d \n", int(_layout[idx].pidx), int(_layout[idx].cnt));
        idx++;
      }
      printf("\n");
    }
    else {
      T idx = 0;
      while(size_t(idx) < _layout.size() && _layout[idx].pidx > -1) {
        printf("%d : %d \n", int(_layout[idx].pidx), int(_layout[idx].cnt));
        idx++;
      }
      printf("\n");
    }
  }

  inline void get_parallel_median_split(const std::vector<S> & vals,
                                        std::vector<bool> & on_left_side)
  {
    int size, rank;
    MPI_Comm_size(_comm, &size);
    MPI_Comm_rank(_comm, &rank);

    std::vector<S> lvals(vals);
    on_left_side.resize(vals.size());

    std::vector<int> perm(vals.size());
    for(size_t i=0; i<perm.size(); i++) perm[i] = i;

    sort_copy(lvals, perm);

    size_t nelem = vals.size();
    double min = nelem ? lvals[0]       : 1e100,
           max = nelem ? lvals[nelem-1] : -1e100;
    MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, _comm);
    MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, _comm);

    double bucket_size = (max*1.05 - min) / double(size);
    std::vector<int> buckets(size_t(size), int(0));
    for(size_t i=0; i<nelem; i++) {
      int idx = (lvals[i] - min) / bucket_size;
      idx = kdpart::clamp(idx, 0, size-1);
      buckets[idx]++;
    }

    std::vector<int> glob_buckets(size_t(size), int(0));
    MPI_Allreduce(buckets.data(), glob_buckets.data(), size, MPI_INT, MPI_SUM, _comm);

    // compute the global number of values and consequently the half number
    int glob_sum = std::accumulate(glob_buckets.begin(), glob_buckets.end(), 0);
    int lhalf = (glob_sum + 1) / 2;
    // figure out the bucket index we have to go through to get the median value
    int dsp = 0, bucket_idx = 0;
    while(bucket_idx < size && (dsp + glob_buckets[bucket_idx]) <= lhalf) {
      dsp += glob_buckets[bucket_idx];
      bucket_idx++;
    }

    // sanity check. we could clamp here, but its better to throw an error since
    // an illegal index should not occur.
    if(bucket_idx < 0 || bucket_idx >= size) {
      fprintf(stderr, "Error: Illigal bucket index !!\n");
    }

    // containers for the values we compute the median on
    std::vector<double> loc_val_bucket, glob_val_bucket;
    // containers for the split decision
    std::vector<short> global_split_bucket, split_bucket(buckets[bucket_idx]);

    loc_val_bucket.resize(buckets[bucket_idx]);
    for(size_t i=0, widx=0; i<nelem; i++) {
      int idx = (lvals[i] - min) / bucket_size;
      idx = kdpart::clamp(idx, 0, size-1);

      if(idx == bucket_idx)
        loc_val_bucket[widx++] = lvals[i];
    }

    std::vector<int> rcnt(size), rdsp(size);
    MPI_Gather(&buckets[bucket_idx], 1, MPI_INT, rcnt.data(), 1, MPI_INT, bucket_idx, _comm);
    kdpart::dsp_from_cnt(rcnt, rdsp);

    if(rank == bucket_idx) {
      glob_val_bucket.resize(glob_buckets[bucket_idx]);
      global_split_bucket.resize(glob_buckets[bucket_idx]);
    }

    // the elements of the split bucket are communicated to the associated rank
    MPI_Gatherv(loc_val_bucket.data(), buckets[bucket_idx], MPI_DOUBLE,
                glob_val_bucket.data(), rcnt.data(), rdsp.data(), MPI_DOUBLE,
                bucket_idx, _comm);


    // the rank owning the bucket where the median split will take place is categorizing
    // his elements into left and right
    if(rank == bucket_idx) {
      size_t bsize = glob_val_bucket.size();

      std::vector<int> bucket_perm(bsize);
      for(size_t i=0; i<bsize; i++) bucket_perm[i] = i;

      sort_copy(glob_val_bucket, bucket_perm);
      int loc_half_idx = lhalf - dsp;

      // sanity check
      if(loc_half_idx < 0 || loc_half_idx >= int(glob_val_bucket.size())) {
        fprintf(stderr, "Error: Illegal local val index!!\n");
      }

      for(int i=0; i <= loc_half_idx; i++)
        global_split_bucket[bucket_perm[i]] = 1;

      for(int i=loc_half_idx+1; i < int(bsize); i++)
        global_split_bucket[bucket_perm[i]] = 0;
    }

    MPI_Scatterv(global_split_bucket.data(), rcnt.data(), rdsp.data(), MPI_SHORT,
                 split_bucket.data(), buckets[bucket_idx], MPI_SHORT,
                 bucket_idx, _comm);

    for(size_t i=0, ridx=0; i<nelem; i++) {
      int pidx = perm[i];
      int idx = (lvals[i] - min) / bucket_size;
      idx = kdpart::clamp(idx, 0, size-1);

      if(idx < bucket_idx)
        on_left_side[pidx] = true;
      else if(idx == bucket_idx)
        on_left_side[pidx] = split_bucket[ridx++] == 1;
      else
        on_left_side[pidx] = false;
    }
  }

  /**
  * @brief Split a parent partition into two child partitions
  *
  * @param parent Parent partition.
  * @param lchild Left child partition.
  * @param rchild Right child partition.
  */
  inline void median_split(kdpart::partition<T,S> parent,
                           kdpart::partition<T,S> & lchild,
                           kdpart::partition<T,S> & rchild)
  {
    kdpart::axis longest_axis = get_longest_axis(parent.box);
    size_t nelem = parent.elems->size();

    std::vector<S> vals(nelem);
    for(size_t i=0; i<nelem; i++) {
      const kdpart::vec3<S> & p = (*parent.elems)[i].ctr;
      S val = S();
      switch(longest_axis) {
        case X: val = p.x; break;
        case Y: val = p.y; break;
        case Z: val = p.z; break;
        case UNSET: break;
      }
      vals[i] = val;
    }

    // we compute the median value in parallel
    std::vector<bool> on_left;
    get_parallel_median_split(vals, on_left);

    // then we split the local elements based on the median
    size_t left_size = 0, right_size = 0;

    for(size_t i=0; i<nelem; i++) {
      if(on_left[i]) left_size++;
      else           right_size++;
    }

    lchild.elems = new std::vector<elem<T,S> >(left_size), rchild.elems = new std::vector<elem<T,S> >(right_size);

    left_size = 0, right_size = 0;
    for(size_t i=0; i<nelem; i++) {
      if(on_left[i]) (*lchild.elems)[left_size++ ] = (*parent.elems)[i];
      else           (*rchild.elems)[right_size++] = (*parent.elems)[i];
    }

    // remove elements of parent partition
    delete parent.elems;
    parent.elems = NULL;

    int lchild_cnt = lchild.elems->size(), rchild_cnt = rchild.elems->size();
    MPI_Allreduce(MPI_IN_PLACE, &lchild_cnt, 1, MPI_INT, MPI_SUM, _comm);
    MPI_Allreduce(MPI_IN_PLACE, &rchild_cnt, 1, MPI_INT, MPI_SUM, _comm);
    lchild.cnt = lchild_cnt, rchild.cnt = rchild_cnt;

    if(lchild_cnt == 0 || rchild_cnt == 0) {
      fprintf(stderr, "Error: Empty partitioning!!\n");
    }

    lchild.box = get_bbox(*lchild.elems);
    rchild.box = get_bbox(*rchild.elems);
  }

  /**
  * @brief Update the paratition layout by splitting a partition
  *
  * @param split_pos   Which partition to split
  */
  inline void update_layout(const T split_pos) {
    assert(size_t(split_pos) < _layout.size());

    T idx_at_split = _layout[split_pos].pidx;
    T start = _curpart - 1, stop = split_pos + 1;

    // increment partition index for partitions after the split
    for(T i = start; i > stop; i--) {
      _layout[i] = _layout[i-1];  // copy partition
      _layout[i].pidx++;          // increment partition index
    }

    // set partition indices for the split
    median_split(_layout[split_pos], _layout[split_pos], _layout[split_pos+1]);
    _layout[split_pos].pidx   = idx_at_split;
    _layout[split_pos+1].pidx = idx_at_split+1;

    // print layout
    // print_layout(split_pos);
  }

  /**
  * @brief Get index of which partition to split next.
  *
  * @return The partition index.
  */
  inline T get_split_pos()
  {
    T idx = 0, max = _layout[0].cnt, maxidx = 0;

    while(idx < _curpart && _layout[idx].cnt > -1) {
      if(max < _layout[idx].cnt) {
        max = _layout[idx].cnt;
        maxidx = idx;
      }
      idx++;
    }

    return maxidx;
  }

  public:
  inline void operator() (const MPI_Comm comm, const std::vector<S> & ctr, const int req_part,
                          std::vector<T> & part_vec)
  {
    _comm = comm;

    assert(ctr.size() % 3 == 0); // coord components need to be a multiple of 3

    // get domain sizes
    long int l_numelem = ctr.size() / 3, g_numelem;
    MPI_Allreduce(&l_numelem, &g_numelem, 1, MPI_LONG, MPI_SUM, _comm);
    assert(g_numelem > 0);

    int size, rank;
    MPI_Comm_size(_comm, &size); MPI_Comm_rank(_comm, &rank);

    T npart = req_part;

    bool redistribute_remainder = false;
    T npart_old = npart;

    /**
    * In the case that the number of processes (thus the number of requested partitions) is
    * not a power of two, we partition with a higher partition number and then evenly
    * distribute the partitions onto the processes. This reduces the difference between
    * partition sizes.
    *
    */
    if(!kdpart::is_power_of_two(npart)) {
      npart = (log(double(npart)) / log(2.0)) + KD_ORDER_INC;
      npart = pow(2, npart);

      // we can always afford to split at least into KD_MIN_SIZE parts. if the initial npart was very small,
      // the computed new npart is still below KD_MIN_SIZE. so we set it explicitely
      if(npart < KD_MIN_SIZE && g_numelem > KD_MIN_SIZE) npart = KD_MIN_SIZE;

      redistribute_remainder = true;
    }

    assert(g_numelem > (long int)npart);

    // initialize the first partition
    _curpart = 1;
    _layout.resize(size_t(npart));
    _layout[0].pidx  = 0;
    _layout[0].elems = new std::vector<kdpart::elem<T,S> >(l_numelem);
    _layout[0].cnt   = g_numelem;

    std::vector<kdpart::elem<T,S> > & elems = *_layout[0].elems;

    // convert the elements into an array of center-points
    for(long int i=0; i<l_numelem; i++) {
      elems[i].ctr.get(ctr.data() + i*3);  // get elem center coord
      elems[i].eidx = i;                   // get elem index
    }

    _layout[0].box = get_bbox(elems);

    /*
     * Main splitting loop: We iterate until we have enough partitions. In each
     * iteration, we split the largest partition into two. Thus we can get to any
     * number of partitions.
     *
     */
    while(_curpart < npart) {
      T split_pos = get_split_pos();
      _curpart++;

      update_layout(split_pos);
    }

    // in case the requested number of partitions was not a power of two,
    // we have computed more partitions than we have processes and we need to re-
    // index the computed partitions into [0, size]
    if(redistribute_remainder) {
      // number of partitions we will at least assign to a process
      T base_size = npart / npart_old;
      // number of partitions we will assign to a subset of processes to reduce
      // the remainder (npart / npart_old)
      T extended_size = base_size + 1;
      T remainder = npart % npart_old;

      // 'remainder' many processes get chunks of size 'extended_size'
      for(T i=0; i<remainder; i++) {
        for(T j=0; j<extended_size; j++)
          _layout[i*extended_size+j].pidx = i;
      }

      // '_layout.size() - remainder' many processes get chunks of size 'base_size'
      for(T i=remainder*extended_size, pidx=remainder; i<T(_layout.size()); i+=base_size, pidx++) {
        for(T j=0; j<base_size; j++)
          _layout[i+j].pidx = pidx;
      }

      // print_layout();
    }

    // assign partition index to the individual elements
    part_vec.assign(l_numelem, T(-1));
    for(const kdpart::partition<T,S> & p : _layout) {
      for(const kdpart::elem<T,S> & e : (*p.elems)) {
        part_vec[e.eidx] = p.pidx;
      }
    }
  }
};
#endif

template<class T, class S>
class sequential_partitioner {

  private:
  // the kdpart members
  std::vector<kdpart::partition<T,S> > _layout;  ///< the partitioning layout
  T _curpart;                              ///< current number of assigned partitions

  /// compute bounding box for a vector of elements
  inline kdpart::bbox<S> get_bbox(const std::vector<elem<T,S> > & elems)
  {
    kdpart::bbox<S> box;
    kdpart::vec3<S> & min = box.bounds[0];
    kdpart::vec3<S> & max = box.bounds[1];
    kdpart::vec3<S> p = elems[0].ctr;

    min.x = p.x, min.y = p.y, min.z = p.z;
    max.x = p.x, max.y = p.y, max.z = p.z;

    for(size_t i=1; i<elems.size(); i++) {
      p = elems[i].ctr;

      if(min.x > p.x) min.x = p.x;
      if(min.y > p.y) min.y = p.y;
      if(min.z > p.z) min.z = p.z;
      if(max.x < p.x) max.x = p.x;
      if(max.y < p.y) max.y = p.y;
      if(max.z < p.z) max.z = p.z;
    }

    return box;
  }

  /// get the longest axis of a bounding box
  inline kdpart::axis get_longest_axis(const kdpart::bbox<S> & box)
  {
    const kdpart::vec3<S> & min = box.bounds[0];
    const kdpart::vec3<S> & max = box.bounds[1];

    S x = fabs(min.x - max.x);
    S y = fabs(min.y - max.y);
    S z = fabs(min.z - max.z);

    return x > y && x > z ? X : y > z ? Y : Z;
  }

  /**
  * @brief Print some details on the layout after a splitting.
  *
  * The location of the last split will be highlighted.
  *
  * @param split_pos The partition where the last splitting occured.
  */
  inline void print_layout(const T split_pos = -1)
  {
    if(split_pos > -1) {
      T idx = 0;
      while(idx < split_pos) {
        printf("%d : %d \n", int(_layout[idx].pidx), int(_layout[idx].cnt));
        idx++;
      }

      printf("----\n");
      printf("%d : %d \n", int(_layout[idx].pidx), int(_layout[idx].cnt));
      printf("%d : %d \n", int(_layout[idx+1].pidx), int(_layout[idx+1].cnt));
      printf("----\n");
      idx += 2;

      while(size_t(idx) < _layout.size() && _layout[idx].pidx > -1) {
        printf("%d : %d \n", int(_layout[idx].pidx), int(_layout[idx].cnt));
        idx++;
      }
      printf("\n");
    }
    else {
      T idx = 0;
      while(size_t(idx) < _layout.size() && _layout[idx].pidx > -1) {
        printf("%d : %d \n", int(_layout[idx].pidx), int(_layout[idx].cnt));
        idx++;
      }
      printf("\n");
    }
  }

  /**
  * @brief Split a parent partition into two child partitions
  *
  * @param parent Parent partition.
  * @param lchild Left child partition.
  * @param rchild Right child partition.
  */
  inline void median_split(kdpart::partition<T,S> parent,
                    kdpart::partition<T,S> & lchild, kdpart::partition<T,S> & rchild)
  {

    kdpart::axis longest_axis = get_longest_axis(parent.box);
    size_t nelem = parent.elems->size();

    // put the (coord, idx) pairs into a vector and sort them
    std::vector<mixed_pair<S,T> > pairs(nelem);
    for(size_t i=0; i<nelem; i++) {
      const kdpart::vec3<S> & p = (*parent.elems)[i].ctr;
      S val = S();
      switch(longest_axis) {
        case X: val = p.x; break;
        case Y: val = p.y; break;
        case Z: val = p.z; break;
        case UNSET: break;
      }

      pairs[i] = {val, T(i)};
    }
    std::sort(pairs.begin(), pairs.end());

    // we then copy the first half into left children and the other half into
    // right children
    size_t lnum  = (nelem + 1) / 2, rnum = nelem - lnum;
    lchild.elems = new std::vector<elem<T,S> >(lnum), rchild.elems = new std::vector<elem<T,S> >(rnum);

    for(size_t i=0; i<lnum; i++) {
      T lidx = pairs[i].v2;
      (*lchild.elems)[i] = (*parent.elems)[lidx];
    }

    for(size_t i=0; i<rnum; i++) {
      T lidx = pairs[lnum + i].v2;
      (*rchild.elems)[i] = (*parent.elems)[lidx];
    }

    // remove elements of parent partition
    delete parent.elems;
    parent.elems = NULL;

    lchild.cnt = lchild.elems->size();
    rchild.cnt = rchild.elems->size();
    lchild.box = get_bbox(*lchild.elems);
    rchild.box = get_bbox(*rchild.elems);
  }

  /**
  * @brief Update the paratition layout by splitting a partition
  *
  * @param split_pos   Which partition to split
  */
  inline void update_layout(const T split_pos) {
    assert(size_t(split_pos) < _layout.size());

    T idx_at_split = _layout[split_pos].pidx;
    T start = _curpart - 1, stop = split_pos + 1;

    // increment partition index for partitions after the split
    for(T i = start; i > stop; i--) {
      _layout[i] = _layout[i-1];  // copy partition
      _layout[i].pidx++;          // increment partition index
    }

    // set partition indices for the split
    median_split(_layout[split_pos], _layout[split_pos], _layout[split_pos+1]);
    _layout[split_pos].pidx   = idx_at_split;
    _layout[split_pos+1].pidx = idx_at_split+1;

    // print layout
    // print_layout(split_pos);
  }

  /**
  * @brief Get index of which partition to split next.
  *
  * @return The partition index.
  */
  inline T get_split_pos()
  {
    T idx = 0, max = _layout[0].cnt, maxidx = 0;

    while(idx < _curpart && _layout[idx].cnt > -1) {
      if(max < _layout[idx].cnt) {
        max = _layout[idx].cnt;
        maxidx = idx;
      }
      idx++;
    }

    return maxidx;
  }

  public:
  inline void operator() (const std::vector<S> & ctr, T npart, std::vector<T> & part_vec)
  {
    assert(ctr.size() % 3 == 0);

    size_t nelem = ctr.size() / 3;
    assert(nelem > size_t(npart));

    bool redistribute_remainder = false;
    T npart_old = npart;

    /**
    * In the case that the number of processes (thus the number of requested partitions) is
    * not a power of two, we partition with a higher partition number and then evenly
    * distribute the partitions onto the processes. This reduces the difference between
    * partition sizes.
    *
    */
    if(!is_power_of_two(npart)) {
      npart = (log(double(npart)) / log(2.0)) + 5.0;
      npart = pow(2, npart);

      // we can always afford to split at least into KD_MIN_SIZE parts. if the initial npart was very small,
      // the computed new npart is still below KD_MIN_SIZE. so we set it explicitely
      if(npart < KD_MIN_SIZE && nelem > KD_MIN_SIZE) npart = KD_MIN_SIZE;

      redistribute_remainder = true;
    }

    // initialize _layout and _cnt
    _layout.resize(size_t(npart));
    part_vec.assign(nelem, T(-1));

    _curpart = 1;
    _layout[0].pidx = 0;
    _layout[0].elems = new std::vector<kdpart::elem<T,S> >(nelem);
    std::vector<kdpart::elem<T,S> > & elems = *_layout[0].elems;

    for(size_t i=0; i<nelem; i++) {
      elems[i].ctr.get(ctr.data() + i*3);  // get elem center coord
      elems[i].eidx = i;                   // get elem index
    }

    _layout[0].box = get_bbox(elems);

    // splitting loop
    while(_curpart < npart) {
      T split_pos = get_split_pos();
      _curpart++;

      update_layout(split_pos);
    }

    if(redistribute_remainder) {
      T base_size = npart / npart_old, extended_size = base_size + 1;
      T remainder = npart % npart_old;

      for(T i=0; i<remainder; i++) {
        for(T j=0; j<extended_size; j++)
          _layout[i*extended_size+j].pidx = i;
      }

      // the rest gets reindexed ascendigly, starting with index "remainder"
      for(T i=remainder*extended_size, pidx=remainder; i<T(_layout.size()); i+=base_size, pidx++) {
        for(T j=0; j<base_size; j++)
          _layout[i+j].pidx = pidx;
      }

      // print_layout();
    }

    // assign partition index to the individual elements
    for(const kdpart::partition<T,S> & p : _layout) {
      for(const kdpart::elem<T,S> & e : (*p.elems)) {
        part_vec[e.eidx] = p.pidx;
      }
    }
  }
};

}

#endif
