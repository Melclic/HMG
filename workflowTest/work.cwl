#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
  w_t:
    type: float
    doc: Doubling time of the cell (min)
    inputBinding:
      position: 1
      prefix: -t
  w_C:
    type: float
    doc: Chromosome replication time (min)
    inputBinding:
      position: 2
      prefix: -c
  w_D:
    type: float
    doc: Cell segregation time (min)
    inputBinding:
      position: 3
      prefix: -d
  w_a:
    type: float
    doc: Age of the cell (0.0<=a<=1.0)
    inputBinding:
      position: 4
      prefix: -a

outputs:
  echoOut:
    type: File
    outputSource: runlibCH/csvOutput

steps:
  runlibCH:
    run: ch.cwl
    in:
      t: w_t
      C: w_C
      D: w_D
      a: w_a
    out: [csvOutput]
