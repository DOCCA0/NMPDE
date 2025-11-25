### 1.1 写出 (1) 的弱形式 (Weak Formulation)

为了推导弱形式，我们将方程乘以测试函数 $v(x)$，并在区间 $\Omega = (0, 1)$ 上进行积分。

原始方程为：
$$\frac{\partial u}{\partial t} - \frac{\partial^2 u}{\partial x^2} + \frac{\partial u}{\partial x} = 0$$
其中 $\kappa=1$。

**步骤：**
1.  **定义空间：**
    * 测试函数空间 $V_0 = H^1_0(0, 1) = \{v \in H^1(0,1) : v(0)=0, v(1)=0\}$。
    * 解空间 $V = \{u \in H^1(0,1) : u(0)=\alpha, u(1)=\beta\}$。
2.  **积分与分部积分（Green公式）：**
    我们将方程两边同乘 $v \in V_0$ 并积分：
    $$\int_0^1 \frac{\partial u}{\partial t} v \, dx + \int_0^1 \left( -\frac{\partial^2 u}{\partial x^2} + \frac{\partial u}{\partial x} \right) v \, dx = 0$$
    对二阶导数项使用分部积分：
    $$\int_0^1 -\frac{\partial^2 u}{\partial x^2} v \, dx = \int_0^1 \frac{\partial u}{\partial x} \frac{\partial v}{\partial x} \, dx - \left[ \frac{\partial u}{\partial x} v \right]_0^1$$
    由于 $v \in V_0$，即 $v(0)=v(1)=0$，边界项 $\left[ \frac{\partial u}{\partial x} v \right]_0^1$ 消失。

**结论（弱形式）：**
寻找 $u(t) \in H^1(0,1)$，满足 $u(0,t)=\alpha, u(1,t)=\beta$，且对于任意测试函数 $v \in H^1_0(0,1)$，成立：
$$\int_0^1 \frac{\partial u}{\partial t} v \, dx + \int_0^1 \frac{\partial u}{\partial x} \frac{\partial v}{\partial x} \, dx + \int_0^1 \frac{\partial u}{\partial x} v \, dx = 0$$
简写为：
$$\left( \frac{\partial u}{\partial t}, v \right) + a(u, v) = 0, \quad \forall v \in H^1_0(0,1)$$
其中双线性形式 $a(u, v) = \int_0^1 (u_x v_x + u_x v) \, dx$。

---

### 1.2 讨论与弱形式相关的双线性形式的强制性 (Coercivity)

我们需要分析双线性形式 $a(u, v)$ 在空间 $H^1_0(0,1)$ 上的强制性（或称椭圆性）。
双线性形式为：
$$a(u, v) = \int_0^1 \left( \frac{\partial u}{\partial x} \frac{\partial v}{\partial x} + \frac{\partial u}{\partial x} v \right) \, dx$$

**分析过程：**
设 $u = v$，其中 $v \in H^1_0(0,1)$。
$$a(v, v) = \int_0^1 \left( \frac{\partial v}{\partial x} \right)^2 \, dx + \int_0^1 \frac{\partial v}{\partial x} v \, dx$$
观察对流项（第二项）：
$$\int_0^1 v' v \, dx = \int_0^1 \frac{1}{2} \frac{d}{dx}(v^2) \, dx = \frac{1}{2} [v(1)^2 - v(0)^2]$$
由于 $v \in H^1_0(0,1)$，即 $v(0)=v(1)=0$，该项为 **0**。
因此：
$$a(v, v) = \int_0^1 |v'|^2 \, dx = |v|_{H^1}^2$$
利用庞加莱不等式（Poincaré inequality），存在常数 $C_p > 0$ 使得对于 $v \in H^1_0(0,1)$ 有 $\|v\|_{L^2} \le C_p |v|_{H^1}$。这意味着半范数 $|v|_{H^1}$ 与全范数 $\|v\|_{H^1}$ 等价。
$$a(v, v) = |v|_{H^1}^2 \ge C \|v\|_{H^1}^2$$
（其中 $C$ 是与庞加莱常数有关的正数）。

**结论：**
双线性形式 $a(\cdot, \cdot)$ 是**强制的 (coercive)**。

---

### 1.3 写出全离散近似格式 (Fully Discrete Approximation)

使用空间上的分段线性有限元（$P1$ 单元）和时间上的隐式向后欧拉法（Backward Euler）。

**定义：**
* **网格：** 将 $[0,1]$ 划分为 $N$ 个单元，节点间距为 $h$。
* **有限元空间：** $V_h \subset H^1(0,1)$ 是由分段线性基函数张成的空间。
    $V_{h,0} = \{v_h \in V_h : v_h(0)=0, v_h(1)=0\}$。
* **时间步：** 时间步长 $\Delta t$，时刻 $t_n = n\Delta t$。$u_h^n$ 表示 $t_n$ 时刻的近似解。

**格式推导：**
时间导数近似为：$\frac{\partial u}{\partial t} \approx \frac{u_h^{n+1} - u_h^n}{\Delta t}$。
在 $t_{n+1}$ 时刻考虑弱形式：

**结论（全离散格式）：**
已知 $u_h^n$，求解 $u_h^{n+1} \in V_h$（满足边界条件 $u_h^{n+1}(0)=\alpha, u_h^{n+1}(1)=\beta$），使得对于任意 $v_h \in V_{h,0}$：
$$\int_0^1 \frac{u_h^{n+1} - u_h^n}{\Delta t} v_h \, dx + a(u_h^{n+1}, v_h) = 0$$
或者写成算子形式：
$$(u_h^{n+1}, v_h) + \Delta t \, a(u_h^{n+1}, v_h) = (u_h^n, v_h), \quad \forall v_h \in V_{h,0}$$

---

### 1.4 当 $\alpha=0, \beta=0$ 时，证明格式的稳定性并给出误差估计

当 $\alpha=\beta=0$ 时，解属于 $V_{h,0}$（齐次Dirichlet边界条件）。

#### 1. 稳定性证明 (Stability Proof)

我们在全离散方程中取测试函数 $v_h = u_h^{n+1}$：
$$\left( \frac{u_h^{n+1} - u_h^n}{\Delta t}, u_h^{n+1} \right) + a(u_h^{n+1}, u_h^{n+1}) = 0$$
由 1.2 的结论可知 $a(u_h^{n+1}, u_h^{n+1}) = |u_h^{n+1}|_{H^1}^2 \ge 0$。
因此：
$$\left( u_h^{n+1} - u_h^n, u_h^{n+1} \right) = -\Delta t \, |u_h^{n+1}|_{H^1}^2 \le 0$$
利用恒等式 $(a-b, a) = \frac{1}{2}\|a\|^2 - \frac{1}{2}\|b\|^2 + \frac{1}{2}\|a-b\|^2$，我们有：
$$\frac{1}{2}\|u_h^{n+1}\|_{L^2}^2 - \frac{1}{2}\|u_h^n\|_{L^2}^2 + \frac{1}{2}\|u_h^{n+1} - u_h^n\|_{L^2}^2 \le 0$$
忽略非负项 $\frac{1}{2}\|u_h^{n+1} - u_h^n\|_{L^2}^2$，得到：
$$\|u_h^{n+1}\|_{L^2}^2 \le \|u_h^n\|_{L^2}^2$$
这意味着数值解的 $L^2$ 范数随时间不增加，即格式是 **无条件稳定** 的（在 $L^2$ 范数意义下）。

#### 2. 误差估计 (Error Estimate)

对于分段线性有限元（空间精度 $O(h^2)$）和向后欧拉法（时间精度 $O(\Delta t)$），标准的先验误差估计为：

**结论：**
$$\| u(t_n) - u_h^n \|_{L^2(\Omega)} \le C (h^2 + \Delta t)$$
其中 $C$ 是与 $h$ 和 $\Delta t$ 无关的常数。