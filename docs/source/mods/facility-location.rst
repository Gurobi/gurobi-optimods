Facility Location
=================

The facility location problem is a fundamental optimization problem in operations
research and logistics. Given a set of customers with known demand and a set of
potential facility locations, the goal is to determine which facilities to open
and how to allocate customer demand to facilities to minimize total costs,
including both fixed costs for opening facilities and variable transportation
costs.

This problem has numerous real-world applications:

- **Warehouse location**: Determining where to build warehouses to serve customers
- **Distribution center planning**: Optimizing the placement of distribution centers
- **Server placement**: Deciding where to locate data centers to serve users
- **Emergency services**: Positioning fire stations, hospitals, or ambulances
- **Retail location**: Choosing store locations to serve a customer base

Problem Specification
---------------------

The capacitated facility location problem can be formally stated as follows:

Given:
  - A set of customers :math:`I`, each with demand :math:`d_i`
  - A set of potential facility locations :math:`J`, each with:

    - Capacity :math:`u_j` (maximum demand it can serve)
    - Fixed cost :math:`f_j` to open the facility

  - Transportation cost :math:`c_{ij}` per unit shipped from facility :math:`j` to customer :math:`i`

Find:
  - Which facilities to open
  - How much to ship from each open facility to each customer

Objective:
  Minimize the total cost (fixed costs + transportation costs)

Subject to:
  - All customer demand must be satisfied
  - Facility capacity limits must be respected
  - Customers can only be served by open facilities
  - Customers may be served by multiple facilities (multiple allocation)

.. dropdown:: Background: Mathematical Model

    The capacitated facility location problem is formulated as a mixed-integer
    linear program (MILP):

    **Decision Variables:**

    - :math:`y_j \in \{0, 1\}`: Binary variable indicating whether facility :math:`j` is opened
    - :math:`x_{ij} \geq 0`: Continuous variable representing the amount shipped from facility :math:`j` to customer :math:`i`

    **Objective Function:**

    .. math::

        \min \sum_{j \in J} f_j y_j + \sum_{i \in I} \sum_{j \in J} c_{ij} x_{ij}

    The objective minimizes the sum of fixed costs for opening facilities and
    transportation costs for serving customers.

    **Constraints:**

    1. **Demand satisfaction:** Each customer's demand must be fully satisfied:

    .. math::

        \sum_{j \in J} x_{ij} = d_i \quad \forall i \in I

    2. **Capacity limits:** The total demand assigned to each facility cannot exceed its capacity:

    .. math::

        \sum_{i \in I} x_{ij} \leq u_j \quad \forall j \in J

    3. **Facility opening:** Customers can only be served by open facilities:

    .. math::

        x_{ij} \leq d_i \cdot y_j \quad \forall i \in I, j \in J

    This is a "big-M" constraint ensuring that if :math:`y_j = 0` (facility closed),
    then :math:`x_{ij} = 0` (no shipments from that facility).

    4. **Variable domains:**

    .. math::

        y_j \in \{0, 1\} \quad \forall j \in J

        x_{ij} \geq 0 \quad \forall i \in I, j \in J

    **Problem Variants:**

    The facility location problem has many variants:

    - **Uncapacitated Facility Location Problem (UFLP):** No capacity limits on facilities
    - **Capacitated Facility Location Problem (CFLP):** Facilities have capacity limits (this implementation)
    - **p-Median Problem:** Open exactly :math:`p` facilities with zero fixed costs
    - **Single-source:** Each customer must be served by exactly one facility
    - **Multi-source:** Customers can be served by multiple facilities (this implementation)

    **References:**

    - Cornuéjols, G., Nemhauser, G. L., & Wolsey, L. A. (1990). The uncapacitated facility location problem.
      In *Discrete location theory* (pp. 119-171). Wiley.
    - Sridharan, R. (1995). The capacitated plant location problem.
      *European Journal of Operational Research*, 87(2), 203-213.

Interface
---------

The ``solve_facility_location`` function accepts three pandas DataFrames as input
describing customers, facilities, and transportation costs. It returns a dictionary
containing the optimal solution.

.. testcode:: facility_location_basic

    import pandas as pd
    from gurobi_optimods.facility_location import solve_facility_location

    # Define customer data: each customer has a demand
    customer_data = pd.DataFrame({
        "customer": ["Boston", "New York", "Philadelphia"],
        "demand": [100.0, 150.0, 120.0]
    })

    # Define facility data: each potential facility has capacity and fixed cost
    facility_data = pd.DataFrame({
        "facility": ["Warehouse A", "Warehouse B", "Warehouse C"],
        "capacity": [200.0, 180.0, 250.0],
        "fixed_cost": [5000.0, 4500.0, 6000.0]
    })

    # Define transportation costs per unit from each facility to each customer
    transportation_cost = pd.DataFrame({
        "customer": ["Boston"] * 3 + ["New York"] * 3 + ["Philadelphia"] * 3,
        "facility": ["Warehouse A", "Warehouse B", "Warehouse C"] * 3,
        "cost": [10.0, 8.0, 15.0,    # Costs to Boston
                 12.0, 7.0, 9.0,      # Costs to New York
                 8.0, 11.0, 6.0]      # Costs to Philadelphia
    })

    # Solve the facility location problem
    result = solve_facility_location(
        customer_data, facility_data, transportation_cost, verbose=False
    )

    # Display which facilities were opened
    opened = result["facilities_opened"]
    print("Facilities opened:")
    for facility in opened[opened > 0.5].index:
        print(f"  - {facility}")

    # Display total cost
    print(f"\nTotal cost: ${result['solution_value']:,.2f}")

.. testoutput:: facility_location_basic
    :options: +NORMALIZE_WHITESPACE

    Facilities opened:
      - Warehouse A
      - Warehouse B

    Total cost: $12,450.00

The solution indicates which facilities to open and how to allocate customer demand:

.. testcode:: facility_location_basic

    # Display customer assignments
    assignments = result["assignments"]
    print("\nCustomer allocations:")
    for _, row in assignments.iterrows():
        print(f"  {row['customer']} <- {row['facility']}: {row['assignment']:.1f} units")

.. testoutput:: facility_location_basic
    :options: +NORMALIZE_WHITESPACE

    Customer allocations:
      Boston <- Warehouse A: 70.0 units
      Boston <- Warehouse B: 30.0 units
      New York <- Warehouse B: 150.0 units
      Philadelphia <- Warehouse A: 120.0 units

Return Value Structure
^^^^^^^^^^^^^^^^^^^^^^

The ``solve_facility_location`` function returns a dictionary with three keys:

- ``solution_value``: The total cost of the optimal solution (float)
- ``facilities_opened``: A pandas Series indexed by facility identifiers, with values of 1.0 for open facilities and 0.0 for closed facilities
- ``assignments``: A pandas DataFrame with columns ``customer``, ``facility``, and ``assignment``, showing the quantity allocated from each facility to each customer (only non-zero allocations are included)

Example: Trade-off Between Fixed and Transportation Costs
----------------------------------------------------------

The facility location problem involves a natural trade-off: opening more facilities
reduces transportation costs (customers are served from nearby facilities) but
increases fixed costs. Let's illustrate this with an example:

.. testcode:: facility_location_tradeoff

    import pandas as pd
    from gurobi_optimods.facility_location import solve_facility_location

    # 4 customers with demand
    customer_data = pd.DataFrame({
        "customer": [1, 2, 3, 4],
        "demand": [20.0, 30.0, 25.0, 35.0]
    })

    # 3 potential facilities with high fixed cost but large capacity
    facility_data = pd.DataFrame({
        "facility": ["A", "B", "C"],
        "capacity": [50.0, 60.0, 70.0],
        "fixed_cost": [100.0, 120.0, 80.0]
    })

    # Transportation costs vary - some facilities are closer to certain customers
    transportation_cost = pd.DataFrame({
        "customer": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
        "facility": ["A", "B", "C"] * 4,
        "cost": [1.0, 5.0, 8.0,     # Customer 1 is close to facility A
                 5.0, 1.0, 7.0,     # Customer 2 is close to facility B
                 8.0, 6.0, 2.0,     # Customer 3 is close to facility C
                 4.0, 2.0, 1.0]     # Customer 4 is close to facility C
    })

    result = solve_facility_location(
        customer_data, facility_data, transportation_cost, verbose=False
    )

    # Analyze the solution
    print(f"Total cost: ${result['solution_value']:.2f}")
    print(f"Facilities opened: {result['facilities_opened'][result['facilities_opened'] > 0.5].index.tolist()}")

    # Break down costs
    opened = result["facilities_opened"]
    fixed_cost_total = (facility_data.set_index("facility")["fixed_cost"] * opened).sum()
    transport_cost_total = result["solution_value"] - fixed_cost_total
    print(f"  Fixed costs: ${fixed_cost_total:.2f}")
    print(f"  Transportation costs: ${transport_cost_total:.2f}")

.. testoutput:: facility_location_tradeoff
    :options: +NORMALIZE_WHITESPACE

    Total cost: $415.00
    Facilities opened: ['B', 'C']
      Fixed costs: $200.00
      Transportation costs: $215.00

In this case, the optimal solution opens only two facilities (B and C), balancing
fixed and transportation costs. Opening all three facilities would cost more
despite potentially reducing some transportation costs.

Example: Capacity Constraints
------------------------------

Facility capacities play a critical role. If a single facility cannot serve all
demand, multiple facilities must be opened:

.. testcode:: facility_location_capacity

    import pandas as pd
    from gurobi_optimods.facility_location import solve_facility_location

    # Two high-demand customers
    customer_data = pd.DataFrame({
        "customer": ["Customer 1", "Customer 2"],
        "demand": [80.0, 90.0]
    })

    # Two facilities with limited capacity
    facility_data = pd.DataFrame({
        "facility": ["Facility 1", "Facility 2"],
        "capacity": [100.0, 100.0],  # Neither can serve all demand alone
        "fixed_cost": [50.0, 60.0]
    })

    # Similar transportation costs
    transportation_cost = pd.DataFrame({
        "customer": ["Customer 1", "Customer 1", "Customer 2", "Customer 2"],
        "facility": ["Facility 1", "Facility 2", "Facility 1", "Facility 2"],
        "cost": [2.0, 2.5, 2.5, 2.0]
    })

    result = solve_facility_location(
        customer_data, facility_data, transportation_cost, verbose=False
    )

    print(f"Facilities opened: {result['facilities_opened'][result['facilities_opened'] > 0.5].index.tolist()}")
    print("\nAllocations:")
    for _, row in result["assignments"].iterrows():
        print(f"  {row['customer']} <- {row['facility']}: {row['assignment']:.0f} units")

.. testoutput:: facility_location_capacity
    :options: +NORMALIZE_WHITESPACE

    Facilities opened: ['Facility 1', 'Facility 2']

    Allocations:
      Customer 1 <- Facility 1: 80 units
      Customer 2 <- Facility 2: 90 units

Both facilities must be opened because the total demand (170 units) exceeds any
single facility's capacity (100 units). In this case, each customer is served
entirely by one facility since the costs are similar.

Example: Fixed Number of Facilities
------------------------------------

In some scenarios, you may want to open exactly a specific number of facilities,
rather than letting the optimizer decide. This is useful when:

- Budget constraints limit the number of facilities that can be built
- Strategic planning requires a specific number of locations
- You want to compare solutions with different numbers of facilities

The ``fixed_facility_count`` parameter enforces this constraint:

.. testcode:: facility_location_fixed_count

    import pandas as pd
    from gurobi_optimods.facility_location import solve_facility_location

    # Define customers and facilities
    customer_data = pd.DataFrame({
        "customer": [1, 2, 3],
        "demand": [40.0, 50.0, 45.0]
    })

    facility_data = pd.DataFrame({
        "facility": ["A", "B", "C"],
        "capacity": [80.0, 90.0, 85.0],
        "fixed_cost": [100.0, 120.0, 110.0]
    })

    transportation_cost = pd.DataFrame({
        "customer": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        "facility": ["A", "B", "C"] * 3,
        "cost": [5.0, 8.0, 12.0,
                 10.0, 4.0, 9.0,
                 12.0, 10.0, 3.0]
    })

    # Solve without constraint (let optimizer decide)
    result_unconstrained = solve_facility_location(
        customer_data, facility_data, transportation_cost, verbose=False
    )

    # Solve with exactly 2 facilities required
    result_constrained = solve_facility_location(
        customer_data, facility_data, transportation_cost,
        fixed_facility_count=2, verbose=False
    )

    # Compare solutions
    print("Unconstrained solution:")
    print(f"  Facilities opened: {result_unconstrained['facilities_opened'].sum():.0f}")
    print(f"  Total cost: ${result_unconstrained['solution_value']:.2f}")

    print("\nConstrained solution (exactly 2 facilities):")
    opened = result_constrained['facilities_opened']
    print(f"  Facilities opened: {opened[opened > 0.5].index.tolist()}")
    print(f"  Total cost: ${result_constrained['solution_value']:.2f}")

.. testoutput:: facility_location_fixed_count
    :options: +NORMALIZE_WHITESPACE

    Unconstrained solution:
      Facilities opened: 3
      Total cost: $865.00

    Constrained solution (exactly 2 facilities):
      Facilities opened: ['B', 'C']
      Total cost: $885.00

In this example, the unconstrained solution opens all 3 facilities for a total
cost of $865, while constraining to exactly 2 facilities yields a slightly
higher cost of $885. This demonstrates how forcing a specific facility count
can lead to suboptimal (but still valid) solutions.

Fairness Constraints
--------------------

The ``max_facilities_per_region`` parameter allows you to limit the concentration
of facilities in specific regions or categories. This is particularly useful
for preventing the unfair distribution of undesirable facilities (such as
landfills, waste treatment plants, or industrial sites) in certain demographic
areas.

**Use Case: Equitable Landfill Distribution**

Consider a scenario where a city needs to locate waste processing facilities.
Without fairness constraints, the optimization might concentrate multiple
facilities in lower-income neighborhoods due to lower land costs. Fairness
constraints ensure a more equitable distribution:

.. testcode:: facility_location_fairness

    import pandas as pd
    from gurobi_optimods.facility_location import solve_facility_location

    # Communities with waste processing needs
    customer_data = pd.DataFrame({
        "customer": ["North", "South", "East", "West"],
        "demand": [45.0, 35.0, 40.0, 30.0]
    })

    # Potential landfill sites in different socioeconomic regions
    facility_data = pd.DataFrame({
        "facility": ["L1", "L2", "M1", "M2", "H1", "H2"],
        "capacity": [60.0, 60.0, 60.0, 60.0, 60.0, 60.0],
        "fixed_cost": [80.0, 85.0, 100.0, 105.0, 120.0, 125.0],
        "region": ["Low-income", "Low-income", "Middle", "Middle", "High-income", "High-income"]
    })

    # Transportation costs ($/ton)
    transportation_cost = pd.DataFrame({
        "customer": ["North"] * 6 + ["South"] * 6 + ["East"] * 6 + ["West"] * 6,
        "facility": ["L1", "L2", "M1", "M2", "H1", "H2"] * 4,
        "cost": [2.0, 2.5, 4.0, 4.5, 6.0, 6.5,
                3.0, 3.5, 2.0, 2.5, 5.0, 5.5,
                2.5, 3.0, 3.0, 3.5, 4.0, 4.5,
                3.0, 3.5, 2.5, 3.0, 4.5, 5.0]
    })

    # Solve without fairness constraints
    result_no_fairness = solve_facility_location(
        customer_data, facility_data, transportation_cost, verbose=False
    )

    # Solve with fairness: at most 1 facility per region
    result_with_fairness = solve_facility_location(
        customer_data, facility_data, transportation_cost,
        max_facilities_per_region=1,
        verbose=False
    )

    # Compare solutions
    print("Without fairness constraints:")
    opened_no_fair = result_no_fairness['facilities_opened']
    for region in ["Low-income", "Middle", "High-income"]:
        region_facilities = facility_data[facility_data['region'] == region]['facility']
        count = opened_no_fair.loc[region_facilities].sum()
        print(f"  {region}: {count:.0f} facilities")
    print(f"  Total cost: ${result_no_fairness['solution_value']:.2f}")

    print("\nWith fairness constraints (max 1 per region):")
    opened_with_fair = result_with_fairness['facilities_opened']
    for region in ["Low-income", "Middle", "High-income"]:
        region_facilities = facility_data[facility_data['region'] == region]['facility']
        count = opened_with_fair.loc[region_facilities].sum()
        opened_list = opened_with_fair.loc[region_facilities]
        opened_list = opened_list[opened_list > 0.5].index.tolist()
        if opened_list:
            print(f"  {region}: {len(opened_list)} facility ({', '.join(opened_list)})")
        else:
            print(f"  {region}: 0 facilities")
    print(f"  Total cost: ${result_with_fairness['solution_value']:.2f}")

.. testoutput:: facility_location_fairness
    :options: +NORMALIZE_WHITESPACE

    Without fairness constraints:
      Low-income: 2 facilities
      Middle: 1 facilities
      High-income: 0 facilities
      Total cost: $617.50

    With fairness constraints (max 1 per region):
      Low-income: 1 facility (L1)
      Middle: 1 facility (M1)
      High-income: 1 facility (H1)
      Total cost: $682.50

In this example, without fairness constraints, the optimizer concentrates
facilities in low-income areas due to lower costs. With fairness constraints
limiting each region to at most 1 facility, the solution distributes facilities
more equitably across all socioeconomic regions. The fairness constraint
increases the total cost by $65 (about 10.5%), which represents the price of
achieving equitable distribution.

**Combining Fairness with Fixed Counts**

Fairness constraints can be combined with ``fixed_facility_count`` to achieve
both an exact number of facilities and equitable distribution:

.. testcode:: facility_location_fairness_combined

    import pandas as pd
    from gurobi_optimods.facility_location import solve_facility_location

    customer_data = pd.DataFrame({
        "customer": [1, 2, 3],
        "demand": [15.0, 15.0, 15.0]
    })

    facility_data = pd.DataFrame({
        "facility": ["A1", "A2", "B1", "B2", "C1", "C2"],
        "capacity": [20.0] * 6,
        "fixed_cost": [10.0, 12.0, 11.0, 13.0, 15.0, 14.0],
        "region": ["R1", "R1", "R2", "R2", "R3", "R3"]
    })

    transportation_cost = pd.DataFrame({
        "customer": [i for i in [1, 2, 3] for _ in range(6)],
        "facility": ["A1", "A2", "B1", "B2", "C1", "C2"] * 3,
        "cost": [1.0] * 18
    })

    result = solve_facility_location(
        customer_data, facility_data, transportation_cost,
        fixed_facility_count=3,
        max_facilities_per_region=1,
        verbose=False
    )

    print(f"Opened exactly {result['facilities_opened'].sum():.0f} facilities")
    print("Distribution:")
    for region in ["R1", "R2", "R3"]:
        region_facilities = facility_data[facility_data['region'] == region]['facility']
        opened = result['facilities_opened'].loc[region_facilities]
        opened_list = opened[opened > 0.5].index.tolist()
        print(f"  {region}: {', '.join(opened_list)}")

.. testoutput:: facility_location_fairness_combined
    :options: +NORMALIZE_WHITESPACE

    Opened exactly 3 facilities
    Distribution:
      R1: A1
      R2: B1
      R3: C2

This ensures exactly 3 facilities are opened with no more than 1 per region,
achieving both the desired count and equitable distribution.

Computing Service Costs from Locations
---------------------------------------

In many real-world applications, you have geographic coordinates for customers
and facilities rather than pre-computed transportation costs. This section shows
how to compute the ``cost`` column in ``transportation_cost`` from coordinate data
using different distance metrics.

**Distance Metrics Overview:**

- **Euclidean (L2 norm)**: Straight-line distance, suitable for air travel or open terrain
- **Manhattan (L1 norm)**: Grid-based distance, suitable for urban street networks
- **Haversine**: Great-circle distance on Earth's surface, for geographic coordinates (lat/lon)

Euclidean Distance Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For facility and customer locations in a planar coordinate system (e.g., within
a city or region), Euclidean distance is often appropriate:

.. testcode:: facility_location_euclidean

    import pandas as pd
    import numpy as np
    from gurobi_optimods.facility_location import solve_facility_location

    # Customer locations (x, y coordinates in km)
    customer_coords = pd.DataFrame({
        "customer": ["C1", "C2", "C3"],
        "x": [10.0, 25.0, 40.0],
        "y": [15.0, 30.0, 10.0],
        "demand": [20.0, 30.0, 25.0]
    })

    # Facility locations (x, y coordinates in km)
    facility_coords = pd.DataFrame({
        "facility": ["F1", "F2"],
        "x": [15.0, 35.0],
        "y": [20.0, 15.0],
        "capacity": [50.0, 50.0],
        "fixed_cost": [100.0, 120.0]
    })

    # Compute Euclidean distances
    transportation_cost = []
    for _, customer in customer_coords.iterrows():
        for _, facility in facility_coords.iterrows():
            distance = np.sqrt(
                (customer['x'] - facility['x'])**2 +
                (customer['y'] - facility['y'])**2
            )
            # Cost is distance times per-km rate (e.g., $2/km)
            cost_per_km = 2.0
            transportation_cost.append({
                "customer": customer['customer'],
                "facility": facility['facility'],
                "cost": distance * cost_per_km
            })

    transportation_cost = pd.DataFrame(transportation_cost)

    # Prepare customer and facility data (without coordinates)
    customer_data = customer_coords[['customer', 'demand']]
    facility_data = facility_coords[['facility', 'capacity', 'fixed_cost']]

    # Solve
    result = solve_facility_location(
        customer_data, facility_data, transportation_cost, verbose=False
    )

    print(f"Total cost: ${result['solution_value']:.2f}")
    print(f"Facilities opened: {result['facilities_opened'][result['facilities_opened'] > 0.5].index.tolist()}")

.. testoutput:: facility_location_euclidean
    :options: +NORMALIZE_WHITESPACE

    Total cost: $1704.92
    Facilities opened: ['F1', 'F2']

Manhattan Distance Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For urban areas with a grid street network, Manhattan distance (L1 norm) better
represents actual travel:

.. testcode:: facility_location_manhattan

    import pandas as pd
    import numpy as np
    from gurobi_optimods.facility_location import solve_facility_location

    # Same locations as before
    customer_coords = pd.DataFrame({
        "customer": ["C1", "C2", "C3"],
        "x": [10.0, 25.0, 40.0],
        "y": [15.0, 30.0, 10.0],
        "demand": [20.0, 30.0, 25.0]
    })

    facility_coords = pd.DataFrame({
        "facility": ["F1", "F2"],
        "x": [15.0, 35.0],
        "y": [20.0, 15.0],
        "capacity": [50.0, 50.0],
        "fixed_cost": [100.0, 120.0]
    })

    # Compute Manhattan distances (L1 norm)
    transportation_cost = []
    for _, customer in customer_coords.iterrows():
        for _, facility in facility_coords.iterrows():
            distance = (
                abs(customer['x'] - facility['x']) +
                abs(customer['y'] - facility['y'])
            )
            cost_per_km = 2.0
            transportation_cost.append({
                "customer": customer['customer'],
                "facility": facility['facility'],
                "cost": distance * cost_per_km
            })

    transportation_cost = pd.DataFrame(transportation_cost)

    customer_data = customer_coords[['customer', 'demand']]
    facility_data = facility_coords[['facility', 'capacity', 'fixed_cost']]

    result = solve_facility_location(
        customer_data, facility_data, transportation_cost, verbose=False
    )

    print(f"Total cost: ${result['solution_value']:.2f}")
    print(f"Facilities opened: {result['facilities_opened'][result['facilities_opened'] > 0.5].index.tolist()}")

.. testoutput:: facility_location_manhattan
    :options: +NORMALIZE_WHITESPACE

    Total cost: $2320.00
    Facilities opened: ['F1', 'F2']

Note that Manhattan distance yields higher transportation costs (2320 vs 1704.92)
since paths must follow the grid rather than straight lines.

Geographic Distance (Haversine) Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For problems spanning large geographic areas, use the Haversine formula to
compute great-circle distances from latitude/longitude coordinates:

.. testcode:: facility_location_haversine

    import pandas as pd
    import numpy as np
    from gurobi_optimods.facility_location import solve_facility_location

    def haversine_distance(lat1, lon1, lat2, lon2):
        """
        Calculate great-circle distance between two points on Earth.
        Returns distance in kilometers.
        """
        # Convert to radians
        lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

        # Haversine formula
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))

        # Earth's radius in kilometers
        r = 6371.0
        return c * r

    # Customer locations (latitude, longitude)
    customer_coords = pd.DataFrame({
        "customer": ["New York", "Chicago", "Houston"],
        "lat": [40.7128, 41.8781, 29.7604],
        "lon": [-74.0060, -87.6298, -95.3698],
        "demand": [100.0, 80.0, 90.0]
    })

    # Facility locations (latitude, longitude)
    facility_coords = pd.DataFrame({
        "facility": ["Memphis", "Dallas"],
        "lat": [35.1495, 32.7767],
        "lon": [-90.0490, -96.7970],
        "capacity": [150.0, 150.0],
        "fixed_cost": [5000.0, 5500.0]
    })

    # Compute Haversine distances
    transportation_cost = []
    for _, customer in customer_coords.iterrows():
        for _, facility in facility_coords.iterrows():
            distance_km = haversine_distance(
                customer['lat'], customer['lon'],
                facility['lat'], facility['lon']
            )
            # Cost per km for long-haul trucking
            cost_per_km = 0.5
            transportation_cost.append({
                "customer": customer['customer'],
                "facility": facility['facility'],
                "cost": distance_km * cost_per_km
            })

    transportation_cost = pd.DataFrame(transportation_cost)

    customer_data = customer_coords[['customer', 'demand']]
    facility_data = facility_coords[['facility', 'capacity', 'fixed_cost']]

    result = solve_facility_location(
        customer_data, facility_data, transportation_cost, verbose=False
    )

    print(f"Total cost: ${result['solution_value']:.2f}")
    opened = result['facilities_opened'][result['facilities_opened'] > 0.5].index.tolist()
    print(f"Facilities opened: {opened}")

.. testoutput:: facility_location_haversine
    :options: +NORMALIZE_WHITESPACE

    Total cost: $142332.82
    Facilities opened: ['Memphis', 'Dallas']

This example shows facility location for distribution centers serving major US
cities, with costs reflecting the actual geographic distances.

**Choosing the Right Metric:**

- Use **Euclidean** for: Small regions, air travel, line-of-sight problems
- Use **Manhattan** for: Urban areas with grid street layouts
- Use **Haversine** for: Large geographic areas, latitude/longitude coordinates
- For real-world accuracy, consider using routing APIs (Google Maps, OSRM) to get actual driving distances

Notes
-----

**Computational Complexity:**
  The facility location problem is NP-hard, meaning there is no known polynomial-time
  algorithm to solve it optimally for all instances. However, modern mixed-integer
  programming solvers like Gurobi can efficiently solve many practical instances.

**Infeasibility:**
  The problem becomes infeasible if the total facility capacity is less than the
  total customer demand. In this case, the function will raise a ``ValueError``.

**Zero Fixed Costs:**
  If all facilities have zero fixed cost, the problem reduces to a transportation
  problem, which is easier to solve.

**Extensions:**
  For more complex scenarios, consider:

  - Adding minimum opening thresholds
  - Including different service levels or time constraints
  - Modeling multi-echelon distribution networks
  - Incorporating uncertainty in demand or costs

See Also
--------

- :doc:`/mods/min-cost-flow`: For pure transportation problems without facility location decisions
- :doc:`/mods/bipartite-matching`: For assignment problems with binary allocations
