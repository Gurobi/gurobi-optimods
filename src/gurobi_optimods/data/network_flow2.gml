graph [
  directed 1
  node [
    id 0
    label "s"
    pos 0
    pos 0
    demand 2
  ]
  node [
    id 1
    label "1"
    pos 1
    pos 0.5
    demand 0
  ]
  node [
    id 2
    label "2"
    pos 1
    pos -0.5
    demand 1
  ]
  node [
    id 3
    label "3"
    pos 2
    pos 0.5
    demand -1
  ]
  node [
    id 4
    label "4"
    pos 2
    pos -0.5
    demand 0
  ]
  node [
    id 5
    label "t"
    pos 3
    pos 0
    demand -2
  ]
  edge [
    source 0
    target 1
    capacity 2
    cost 9
  ]
  edge [
    source 0
    target 2
    capacity 2
    cost 7
  ]
  edge [
    source 1
    target 3
    capacity 1
    cost 1
  ]
  edge [
    source 2
    target 3
    capacity 1
    cost 10
  ]
  edge [
    source 2
    target 4
    capacity 2
    cost 6
  ]
  edge [
    source 3
    target 5
    capacity 2
    cost 1
  ]
  edge [
    source 4
    target 5
    capacity 2
    cost 1
  ]
]
