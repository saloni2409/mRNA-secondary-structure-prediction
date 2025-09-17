"""
Module for building the QUBO formulation for RNA secondary structure prediction.

This module translates the RNA folding problem into a Quadratic Unconstrained
Binary Optimization (QUBO) problem, which can then be solved by quantum
or classical optimizers.

--- Quadratic Unconstrained Binary Optimization (QUBO) ---

In the context of RNA folding, the problem is modeled by assigning a binary
variable to each potential structural component. Here, the component is a
"quartet," a stable stack of two adjacent base pairs.

- Binary Variable: For each possible quartet `q_i`, we define a binary variable.
  - `q_i = 1` means "we are choosing to include this specific quartet in our
    final folded structure."
  - `q_i = 0` means we are not.

- Objective Function: The goal is to find the combination of 1s and 0s (the
  selection of quartets) that results in the minimum possible total free energy
  of the folded structure. This is known as finding the Minimum Free Energy (MFE)
  structure. The objective function is a mathematical formula that calculates this
  total energy for any given set of choices.

  Total Energy = (Sum of individual quartet energies) + (Sum of interaction rewards/penalties)

--- From Constrained to Unconstrained Problem ---

A key rule in RNA secondary structure is that base pairs cannot "cross." In a
formal constrained optimization, this would be a separate rule like `q_i + q_j <= 1`
for any two crossing quartets `i` and `j`.

To make this an *unconstrained* problem suitable for QUBO, we move this rule
into the objective function itself. We add a large positive energy penalty (`t`)
to the objective function for any two incompatible (crossing) quartets. The
term `t * q_i * q_j` is added to the total energy. If an optimizer tries to
select both `q_i` and `q_j` (so `q_i=1` and `q_j=1`), the total energy gets a
massive penalty, making that solution highly unfavorable. The optimizer will
naturally discard it in favor of lower-energy, valid structures.

--- QUBO Objective Function Components ---

The QUBO objective function is a summation of linear and quadratic terms.

1.  **Linear Terms (Individual Costs/Rewards):**
    Each potential quartet has an intrinsic stability, corresponding to a certain
    amount of free energy, `e_qi`. These values are typically negative (favorable)
    and are determined experimentally. The linear part of the objective function
    is the sum of the energies of all chosen quartets.

    `F_linear = Σ (e_qi * q_i)` for all possible quartets `q_i`

2.  **Quadratic Terms (Interaction Costs/Rewards):**
    These terms represent the bonus or penalty that arises from the interaction
    between two choices (i.e., selecting two different quartets).

    -   **Stacking Rewards:** When two quartets are stacked directly on top of each
        other, they create a more stable helix. This is represented by an energy
        reward `r` (a negative value). The term `r * q_i * q_j` is added for each
        pair of stackable quartets `i` and `j`. This reward is only applied if
        both quartets are chosen (`q_i=1` and `q_j=1`).

        `F_reward = Σ (r * q_i * q_j)` for all stackable pairs `(q_i, q_j)`

    -   **Crossing Penalties:** As described above, a large penalty `t` (a positive
        value) is added for any pair of quartets that cross each other.

        `F_penalty = Σ (t * q_i * q_j)` for all crossing pairs `(q_i, q_j)`

--- Final QUBO Formulation Implemented ---

The complete objective function to be minimized is:

`Min F(q) = F_linear + F_reward + F_penalty`
`Min F(q) = Σ(e_qi * q_i) + Σ(r * q_i * q_j) + Σ(t * q_i * q_j)`

This module's `build_rna_qubo` function constructs this exact objective
function as a `qiskit_optimization.problems.QuadraticProgram` object.

Note: This model can be extended with other penalties, such as for a helix
ending in a (U,A) pair, by adding more linear and quadratic terms to the
objective function.
"""
# Import the necessary Qiskit class
from qiskit_optimization.problems import QuadraticProgram

# --- 1. Helper Functions for Pre-processing ---

def is_valid_pair(b1, b2):
    """Checks if two bases form a valid pair (Watson-Crick or wobble)."""
    return tuple(sorted((b1, b2))) in [('A', 'U'), ('C', 'G'), ('G', 'U')]

def generate_quartets(sequence):
    """Generates all possible valid quartets from an RNA sequence."""
    n = len(sequence)
    quartets = []
    # A quartet is two stacked pairs (i, j) and (i+1, j-1)
    # The variable is represented as a tuple: (i, j, i+1, j-1)
    for i in range(n):
        for j in range(i + 4, n): # Ensure a minimal loop size
            # Check if the outer pair (i, j) is valid
            if is_valid_pair(sequence[i], sequence[j]):
                # Check if the inner pair (i+1, j-1) is valid
                if is_valid_pair(sequence[i+1], sequence[j-1]):
                    quartets.append((i, j, i+1, j-1))
    return quartets

def are_quartets_crossing(q1, q2):
    """Checks if two quartets contain crossing pairs."""
    # q1 = (i, j, i+1, j-1) and q2 = (k, l, k+1, l-1)
    i, j = q1[0], q1[1]
    k, l = q2[0], q2[1]
    # Simple crossing condition: i < k < j < l or k < i < l < j
    return (i < k and k < j and j < l) or \
           (k < i and i < l and l < j)

def are_quartets_stacked(q1, q2):
    """Checks if q2 is stacked directly on top of q1."""
    # q1 = (i, j, i+1, j-1) and q2 should be (i+2, j-2, ...)
    return q1[2] + 1 == q2[0] and q1[3] - 1 == q2[1]


# --- 2. The Main QUBO Building Function ---

def build_rna_qubo(sequence, energy_data, stacking_reward, crossing_penalty):
    """
    Builds the QUBO for RNA folding as a Qiskit QuadraticProgram.

    Args:
        sequence (str): The RNA sequence (e.g., "GCGAUAGCGC").
        energy_data (dict): A dictionary mapping quartet tuples to their free energy.
        stacking_reward (float): The energy reward for stacked quartets (should be negative).
        crossing_penalty (float): The energy penalty for crossing quartets (should be large and positive).

    Returns:
        QuadraticProgram: The Qiskit object representing the QUBO.
    """
    # This pre-processing step generates the variables of our problem [cite: 152]
    all_quartets = generate_quartets(sequence)
    
    # Initialize the Quadratic Program for our QUBO
    qp = QuadraticProgram("RNA-Folding-QUBO")

    # --- 3. Add Binary Variables ---
    # For each possible quartet, create a binary variable q_i 
    for quartet in all_quartets:
        qp.binary_var(name=str(quartet))

    # --- 4. Build the Objective Function ---
    linear_terms = {}
    quadratic_terms = {}

    # Add Linear Terms: The base energy for each chosen quartet 
    # This corresponds to the sum: Σ e_qi * q_i
    for quartet in all_quartets:
        q_name = str(quartet)
        linear_terms[q_name] = energy_data.get(quartet, 0)

    # Add Quadratic Terms: Rewards and penalties for interactions
    for i in range(len(all_quartets)):
        for j in range(i + 1, len(all_quartets)):
            q1 = all_quartets[i]
            q2 = all_quartets[j]
            q1_name = str(q1)
            q2_name = str(q2)

            # Add stacking reward: r * q_i * q_j 
            if are_quartets_stacked(q1, q2) or are_quartets_stacked(q2, q1):
                quadratic_terms[(q1_name, q2_name)] = stacking_reward
            
            # Add crossing penalty: t * q_i * q_j [cite: 153, 171]
            if are_quartets_crossing(q1, q2):
                # We add to the existing value in case a pair has multiple interactions
                current_val = quadratic_terms.get((q1_name, q2_name), 0)
                quadratic_terms[(q1_name, q2_name)] = current_val + crossing_penalty
    
    # Set the objective in the QuadraticProgram object
    qp.minimize(linear=linear_terms, quadratic=quadratic_terms)
    
    return qp


# --- 5. Example Usage ---
if __name__ == "__main__":
    # A simple RNA sequence
    rna_sequence = "GCGAUAGCGC"
    
    # Dummy data for the example (in a real scenario, this comes from a database)
    # Let's say we found two possible quartets:
    q1_tuple = (0, 9, 1, 8)  # Pairs G(0)-C(9) and C(1)-G(8)
    q2_tuple = (2, 7, 3, 6)  # Pairs G(2)-C(7) and A(3)-G(6) - Let's assume A-G is valid for this example
    
    dummy_energy_data = {
        q1_tuple: -3.3,
        q2_tuple: -2.1,
    }
    
    # Define reward and penalty values
    r = -1.5  # Stacking reward
    t = 100.0   # Crossing penalty
    
    # Build the QUBO
    rna_qubo = build_rna_qubo(
        sequence=rna_sequence,
        energy_data=dummy_energy_data,
        stacking_reward=r,
        crossing_penalty=t
    )
    
    # Print the QUBO in a human-readable format
    print("--- Generated QUBO ---")
    print(rna_qubo.export_as_lp_string())