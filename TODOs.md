## General

* Use `read_surf_info` where surface data is read.
* Use `apply_normal_orientation` where surface normal vector direction is important.
* change how co-located vertices are found from the current hashmap implementation
  to a kdtree-based implementation.

## Resampling

* Check collapsing criteria from literature.
* Concept on how to do approach uniform refinement in a mesh region:
    1) Compute interface between region to refine and its complement.
    2) Refine interface edges on complement.
    3) Do uniform mesh refinement on region.
    4) Merge region and complement.
