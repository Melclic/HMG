#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['']
inputs:
  t:
    type: float
    doc: Doubling time of the cell (min)
    inputBinding:
      position: 1
      prefix: -t
  C:
    type: float
    doc: Chromosome replication time (min)
    inputBinding:
      position: 2
      prefix: -c
  D:
    type: float
    doc: Cell segregation time (min)  
    inputBinding:
      position: 3
      prefix: -d
  a:
    type: float
    doc: Age of the cell (0.0<=a<=1.0) 
    inputBinding:
      position: 4
      prefix: -a
outputs:
  csvOutput:
    type: File
    outputBinding:
      glob: result.csv
    doc: Output csv file
