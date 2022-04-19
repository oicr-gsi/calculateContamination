digraph calculateContamination {

  subgraph cluster_0 {
    shape=Mrecord
    style=filled;
    color=lightgrey;
    node [style=filled,color=white];
    tumorOnlyMetrics;
    getMetrics
  }

  subgraph cluster_1 {
    node [style=filled];
    b1;
    color=blue
  }
  
  subgraph cluster_2 {
    node [style=filled];
    b3;
    color=blue
  }  
  inputGroups -> bwaMem;
  bwaMem -> bams
  bams -> tumorOnlyMetrics [headlabel="If array length is 1" labeldistance=5 labelangle=120]
  bams -> getMetrics [headlabel="If array length is 2" labeldistance=5 labelangle=-120]
  tumorOnlyMetrics -> b1;
  getMetrics -> b3;
  b1 -> end;
  b3 -> end;

  bams [shape=diamond, label="Bam files", style=rounded, color=orange];
  inputGroups [shape=diamond,label="FastQ files", style=rounded, color=orange];
  end [shape=record, label="Table of contamination metrics"];
  b1 [label="Run non-matched analysis", shape=Mrecord]
  b3 [label="Run tumor/normal analysis", shape=Mrecord]
}