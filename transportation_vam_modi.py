"""
Transportation Optimization using:
1. Vogel's Approximation Method (VAM)
2. MODI (Modified Distribution) Method

The program accepts user input and solves any BALANCED transportation
problem with positive supply and demand values.

Workflow:
Input -> VAM -> Initial Basic Feasible Solution -> MODI ->
Optimality Test -> Improvement -> Optimal Solution
"""

from collections import deque


def read_vector(prompt, size):
    """Read exactly 'size' non-negative integer values."""
    while True:
        try:
            values = list(map(int, input(prompt).split()))

            if len(values) != size:
                print(f"Please enter exactly {size} values.")
                continue

            if any(x <= 0 for x in values):
                print("All supply and demand values must be positive.")
                continue

            return values

        except ValueError:
            print("Please enter integers separated by spaces.")


def print_matrix(title, matrix):
    print(f"\n{title}")
    for row in matrix:
        print(" ".join(f"{x:6}" for x in row))


def calculate_cost(cost, allocation):
    return sum(
        cost[i][j] * allocation[i][j]
        for i in range(len(cost))
        for j in range(len(cost[0]))
    )


# ------------------------------------------------------------------
# VAM
# ------------------------------------------------------------------

def vogel_approximation(cost, supply, demand):
    """Find an Initial Basic Feasible Solution using VAM."""

    m, n = len(supply), len(demand)

    s = supply[:]
    d = demand[:]

    allocation = [[0] * n for _ in range(m)]

    while sum(s) > 0:

        candidates = []

        # Row penalties
        for i in range(m):

            if s[i] == 0:
                continue

            values = sorted(
                cost[i][j]
                for j in range(n)
                if d[j] > 0
            )

            penalty = (
                values[1] - values[0]
                if len(values) >= 2
                else values[0]
            )

            candidates.append(
                (penalty, values[0], 1, -i, i, 'row')
            )

        # Column penalties
        for j in range(n):

            if d[j] == 0:
                continue

            values = sorted(
                cost[i][j]
                for i in range(m)
                if s[i] > 0
            )

            penalty = (
                values[1] - values[0]
                if len(values) >= 2
                else values[0]
            )

            candidates.append(
                (penalty, values[0], 0, -j, j, 'col')
            )

        # Largest penalty.
        # Ties: lower minimum cost, row before column,
        # then smaller index.
        _, _, _, _, index, kind = max(candidates)

        if kind == 'row':

            i = index

            j = min(
                (j for j in range(n) if d[j] > 0),
                key=lambda j: (cost[i][j], j)
            )

        else:

            j = index

            i = min(
                (i for i in range(m) if s[i] > 0),
                key=lambda i: (cost[i][j], i)
            )

        quantity = min(s[i], d[j])

        allocation[i][j] = quantity

        s[i] -= quantity
        d[j] -= quantity

    return allocation


# ------------------------------------------------------------------
# Degeneracy / Basis
# ------------------------------------------------------------------

def create_basis(allocation, cost):
    """
    A transportation basis must contain m+n-1 cells.

    Positive allocations are first included. If the solution is
    degenerate, zero-allocation cells are added without creating
    a cycle.
    """

    m, n = len(allocation), len(allocation[0])

    basis = {
        (i, j)
        for i in range(m)
        for j in range(n)
        if allocation[i][j] > 0
    }

    parent = list(range(m + n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)

        if a == b:
            return False

        parent[a] = b
        return True

    # Positive allocation cells must form a forest.
    for i, j in sorted(basis):
        if not union(i, m + j):
            raise ValueError("Invalid transportation basis.")

    # Add zero cells only when they do not create a cycle.
    candidates = sorted(
        (cost[i][j], i, j)
        for i in range(m)
        for j in range(n)
        if (i, j) not in basis
    )

    for _, i, j in candidates:

        if len(basis) == m + n - 1:
            break

        if union(i, m + j):
            basis.add((i, j))

    if len(basis) != m + n - 1:
        raise ValueError("Could not construct a valid basis.")

    return basis


# ------------------------------------------------------------------
# MODI Potentials
# ------------------------------------------------------------------

def calculate_potentials(cost, basis):
    """
    For every basic cell:
        u[i] + v[j] = cost[i][j]

    We set u[0] = 0 and calculate all other potentials.
    """

    m, n = len(cost), len(cost[0])

    u = [None] * m
    v = [None] * n

    u[0] = 0

    queue = deque([("row", 0)])

    # Bipartite graph:
    # row nodes <-> column nodes
    graph = [[] for _ in range(m + n)]

    for i, j in basis:

        row_node = i
        col_node = m + j

        graph[row_node].append((col_node, cost[i][j]))
        graph[col_node].append((row_node, cost[i][j]))

    visited = [False] * (m + n)
    visited[0] = True

    while queue:

        node_type, node = queue.popleft()

        for next_node, cell_cost in graph[node]:

            if visited[next_node]:
                continue

            visited[next_node] = True

            if next_node < m:
                # Moving from column to row:
                # u[i] = c[i][j] - v[j]
                j = node - m
                u[next_node] = cell_cost - v[j]

                queue.append(("row", next_node))

            else:
                # Moving from row to column:
                # v[j] = c[i][j] - u[i]
                i = node
                v[next_node - m] = cell_cost - u[i]

                queue.append(("col", next_node))

    if any(x is None for x in u) or any(x is None for x in v):
        raise ValueError("Basis is disconnected; MODI cannot continue.")

    return u, v


# ------------------------------------------------------------------
# MODI Loop
# ------------------------------------------------------------------

def find_cycle(m, n, basis, entering):
    """
    Find the unique closed loop created when 'entering' is added
    to a tree basis.

    The loop alternates between rows and columns.
    """

    start = entering[0]
    target = m + entering[1]

    graph = [[] for _ in range(m + n)]

    for i, j in basis:

        a = i
        b = m + j

        graph[a].append((b, (i, j)))
        graph[b].append((a, (i, j)))

    # Find path from entering row to entering column.
    previous = [None] * (m + n)
    previous[start] = (-1, None)

    queue = deque([start])

    while queue:

        node = queue.popleft()

        if node == target:
            break

        for next_node, cell in graph[node]:

            if previous[next_node] is None:

                previous[next_node] = (node, cell)
                queue.append(next_node)

    if previous[target] is None:
        raise ValueError("Could not find an improvement loop.")

    path = []

    node = target

    while node != start:

        prev_node, cell = previous[node]

        path.append(cell)
        node = prev_node

    path.reverse()

    # Entering cell is '+'.
    # Remaining cells are the path in reverse direction.
    return [entering] + list(reversed(path))


# ------------------------------------------------------------------
# MODI Optimization
# ------------------------------------------------------------------

def modi_optimize(cost, allocation):
    """Improve the VAM solution until it is optimal."""

    m, n = len(cost), len(cost[0])

    basis = create_basis(allocation, cost)

    iteration = 1

    while True:

        u, v = calculate_potentials(cost, basis)

        # Calculate opportunity costs:
        #
        # delta = cost - (u + v)
        #
        # For minimization:
        # all delta >= 0 => optimal.
        deltas = []

        for i in range(m):

            for j in range(n):

                if (i, j) not in basis:

                    delta = cost[i][j] - u[i] - v[j]

                    deltas.append((delta, i, j))

        # Most negative delta enters the basis.
        delta, er, ec = min(deltas)

        print(f"\nMODI Iteration {iteration}")
        print("U potentials:", u)
        print("V potentials:", v)

        print("Opportunity costs:")
        for dlt, r, c in deltas:
            print(f"  Delta({r + 1},{c + 1}) = {dlt}")

        # For a minimization problem, all delta >= 0
        # means the current solution is optimal.
        if delta >= 0:

            print("All opportunity costs are >= 0.")
            print("Optimal solution reached.")

            return allocation

        print(f"Entering cell = ({er + 1}, {ec + 1})")
        print(f"Negative opportunity cost = {delta}")

        # Find closed loop.
        cycle = find_cycle(m, n, basis, (er, ec))

        print("Improvement loop:")

        for k, (r, c) in enumerate(cycle):
            sign = "+" if k % 2 == 0 else "-"
            print(f"  {sign} ({r + 1}, {c + 1})")

        plus_cells = cycle[0::2]
        minus_cells = cycle[1::2]

        # theta = smallest allocation on minus cells.
        theta = min(
            allocation[r][c]
            for r, c in minus_cells
        )

        print(f"Theta = {theta}")

        # Add theta to all plus cells.
        for r, c in plus_cells:

            allocation[r][c] += theta

            if (r, c) not in basis:
                basis.add((r, c))

        # Subtract theta from ALL minus cells.
        zero_cells = []

        for r, c in minus_cells:

            allocation[r][c] -= theta

            if allocation[r][c] == 0:
                zero_cells.append((r, c))

        # One zero cell leaves the basis.
        # If multiple cells become zero, the others remain
        # as zero basic cells to avoid changing the basis size.
        if zero_cells:
            basis.remove(zero_cells[0])

        iteration += 1

        if iteration > 1000:
            raise ValueError("Too many MODI iterations.")


# ------------------------------------------------------------------
# Main Program
# ------------------------------------------------------------------

def solve_transportation():

    print("=" * 65)
    print("TRANSPORTATION OPTIMIZATION USING VAM + MODI")
    print("=" * 65)

    # ---------------- INPUT ----------------

    R = int(input("\nEnter number of sources: "))
    C = int(input("Enter number of destinations: "))

    if R < 2 or C < 2:
        print("There must be at least 2 sources and 2 destinations.")
        return

    print("\nEnter transportation cost matrix:")

    cost = []

    for i in range(R):

        while True:

            try:
                row = list(
                    map(
                        int,
                        input(f"Source {i + 1}: ").split()
                    )
                )

                if len(row) != C:
                    print(f"Enter exactly {C} costs.")
                    continue

                if any(x < 0 for x in row):
                    print("Transportation costs cannot be negative.")
                    continue

                cost.append(row)
                break

            except ValueError:
                print("Please enter integers only.")

    supply = read_vector(
        f"\nEnter supply of {R} sources: ",
        R
    )

    demand = read_vector(
        f"Enter demand of {C} destinations: ",
        C
    )

    # Balanced transportation problem.
    if sum(supply) != sum(demand):

        print("\nERROR")
        print("Total supply and total demand must be equal.")
        print("Supply =", sum(supply))
        print("Demand =", sum(demand))
        print("Please enter a balanced transportation problem.")

        return

    print("\nProblem is balanced.")
    print("Total supply =", sum(supply))
    print("Total demand =", sum(demand))

    # ---------------- VAM ----------------

    print("\n" + "=" * 65)
    print("1. VOGEL'S APPROXIMATION METHOD (VAM)")
    print("=" * 65)

    allocation = vogel_approximation(
        cost,
        supply,
        demand
    )

    print_matrix(
        "Initial Basic Feasible Solution (VAM):",
        allocation
    )

    print(
        "\nVAM transportation cost =",
        calculate_cost(cost, allocation)
    )

    # ---------------- MODI ----------------

    print("\n" + "=" * 65)
    print("2. MODI OPTIMALITY TEST AND IMPROVEMENT")
    print("=" * 65)

    allocation = modi_optimize(
        cost,
        allocation
    )

    # ---------------- FINAL ----------------

    print("\n" + "=" * 65)
    print("3. FINAL OPTIMAL SOLUTION")
    print("=" * 65)

    print_matrix(
        "Optimal Shipment Allocation:",
        allocation
    )

    minimum_cost = calculate_cost(
        cost,
        allocation
    )

    print("\nMinimum Transportation Cost =", minimum_cost)

    print("\nShipment Plan:")

    for i in range(R):

        for j in range(C):

            if allocation[i][j] > 0:

                print(
                    f"Source {i + 1} -> Destination {j + 1}"
                    f" : {allocation[i][j]} units"
                )

    # ---------------- Verification ----------------

    print("\n" + "=" * 65)
    print("4. FEASIBILITY CHECK")
    print("=" * 65)

    valid = True

    for i in range(R):

        row_total = sum(allocation[i])

        print(
            f"Source {i + 1}: "
            f"allocated = {row_total}, "
            f"supply = {supply[i]}"
        )

        if row_total != supply[i]:
            valid = False

    for j in range(C):

        col_total = sum(
            allocation[i][j]
            for i in range(R)
        )

        print(
            f"Destination {j + 1}: "
            f"received = {col_total}, "
            f"demand = {demand[j]}"
        )

        if col_total != demand[j]:
            valid = False

    if valid:
        print("\nAll supply and demand constraints are satisfied.")
    else:
        print("\nERROR: Final allocation is not feasible.")


if __name__ == "__main__":
    solve_transportation()
