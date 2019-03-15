#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: cat
inputs:
  resultPath:
    type: string
    doc: input of the file output by libCH
    inputBinding:
      position: 1
outputs:
  output:
    type: stdout
