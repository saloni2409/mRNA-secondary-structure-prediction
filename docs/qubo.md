I'm sorry, I cannot create a file that is directly downloadable to a specific platform like Canvas. However, I can provide the content in Markdown format, which you can easily copy and paste into a file on your computer with a `.md` extension. This file can then be uploaded to Canvas.

---

### **Quadratic Unconstrained Binary Optimization (QUBO) for RNA Folding**

**Quadratic Unconstrained Binary Optimization (QUBO)** is a method used to solve optimization problems by finding the best possible combination of binary variables—variables that can only be 0 or 1. In the context of RNA folding, each potential **quartet** (a stable stack of two base pairs) is a binary variable. A value of $q_i = 1$ means that a specific quartet is included in the final folded structure, while $q_i = 0$ means it is not.

The goal is to find the combination of quartets that results in the **Minimum Free Energy (MFE)** structure. This is done by minimizing an **objective function**, a mathematical formula that calculates a score, or total free energy, for any given set of choices.

---

### **Handling Constraints with Penalties**

In many optimization problems, certain choices are incompatible. For example, in RNA folding, some quartets are physically incompatible because they might use the same base or cross each other in an invalid way. In a standard optimization problem, these would be handled as separate rules or **constraints**.

In QUBO, these constraints are incorporated directly into the objective function as a large **energy penalty**. By adding a significant penalty for incompatible choices, the optimizer will automatically avoid any solution that contains these invalid structures, making them "too expensive". For instance, a rule that two crossing quartets, $q_i$ and $q_j$, cannot both be selected is handled by adding a large penalty, $t$, to the objective function if both $q_i = 1$ and $q_j = 1$.

---

### **Components of the QUBO Objective Function**

The QUBO objective function is a sum of **linear terms** and **quadratic terms**.

* **Linear Terms (Individual Costs/Rewards):** These terms represent the intrinsic energy of each individual quartet. Every quartet has a certain amount of free energy, determined experimentally and available in databases. The sum of these individual energies for all chosen quartets is calculated using the formula:
    $$F_{linear} = \sum_{q_i \in Q} e_{q_i} q_i$$
    Here, $Q$ is the set of all possible valid quartets, $e_{q_i}$ is the energy of quartet $q_i$, and the term is only included if $q_i=1$.

* **Quadratic Terms (Interaction Costs/Rewards):** These terms account for the bonus or penalty that arises from the interaction between two choices.

    * **Stacking Rewards:** When two quartets are stacked consecutively, there is an energy reward (a negative energy value) that lowers the total energy. This is calculated as:
        $$F_{reward} = r \sum_{q_i \in Q, q_j \in QS(q_i)} q_i q_j$$
        Here, $r$ is the reward value, and the term $q_i q_j$ is only 1 if both quartets $q_i$ and $q_j$ are chosen. $QS(q_i)$ represents all quartets that can be stacked with quartet $q_i$.

    * **Penalties:** Penalties are added for unfavorable structures, such as a helix ending in a (U,A) pair. A penalty value, $p$, is applied if a quartet $q_i$ is chosen ($q_i = 1$) and it is at the end of a stack, meaning the next potential stacking quartet $q_j$ is not chosen ($q_j = 0$). This is represented by the formula:
        $$F_{UA\_penalty} = p \sum_{q_i \in Q, q_j \in Q_{UA}} q_i(1 - q_j)$$
        $Q_{UA}$ is the set of quartets with a (U,A) end pair.

    * **Crossing Penalties:** A large penalty, $t$, is added to the objective function if two quartets, $q_i$ and $q_j$, that cannot both be selected are chosen. This ensures the optimizer avoids these invalid combinations. This is calculated as:
        $$F_{constraint\_penalty} = t \sum_{(q_i, q_j) \in Q_C} q_i q_j$$
        $Q_C$ is the set of all pairs of quartets that cross each other.

---

### **The Full Objective Function**

The total objective function to be minimized is the sum of all these terms:

$Min \ F(q) = F_{linear} + F_{reward} + F_{UA\_penalty} + F_{constraint\_penalty}$

$Min \ F(q) = \sum_{q_i \in Q} e_{q_i} q_i + r \sum_{q_i, q_j} q_i q_j + p \sum_{q_i, q_j} q_i(1 - q_j) + t \sum_{(q_i, q_j) \in Q_C} q_i q_j$