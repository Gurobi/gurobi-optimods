graph [
  directed 1
  node [
    id 0
    label "0"
    demand 20
  ]
  node [
    id 1
    label "1"
    demand 0
  ]
  node [
    id 2
    label "2"
    demand 0
  ]
  node [
    id 3
    label "3"
    demand -5
  ]
  node [
    id 4
    label "4"
    demand -15
  ]
  edge [
    source 0
    target 1
    capacity 15
    cost 4
  ]
  edge [
    source 0
    target 2
    capacity 8
    cost 4
  ]
  edge [
    source 1
    target 3
    capacity 4
    cost 2
  ]
  edge [
    source 1
    target 2
    capacity 20
    cost 2
  ]
  edge [
    source 1
    target 4
    capacity 10
    cost 6
  ]
  edge [
    source 2
    target 3
    capacity 15
    cost 1
  ]
  edge [
    source 2
    target 4
    capacity 5
    cost 3
  ]
  edge [
    source 3
    target 4
    capacity 20
    cost 2
  ]
  edge [
    source 4
    target 2
    capacity 4
    cost 3
  ]
]
