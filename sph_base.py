import taichi as ti
import numpy as np


@ti.data_oriented
class SPHBase:
    def __init__(self, particle_system):
        self.ps            = particle_system
        self.g             = -1.80  # Gravity
        self.viscosity     = 0.08  # viscosity
        self.density_0     = 1000.0  # reference density
        self.mass          = self.ps.m_V * self.density_0
        self.dt            = ti.field(ti.f32, shape=())
        self.dt[None]      = 100e-4

        self.up_area       = ti.field(ti.f32, shape=(2,))
        self.up_area[0]    = 4-2/3
        self.up_area[1]    = 4+2/3

    @ti.func
    def cubic_kernel(self, r_norm):
        res = ti.cast(0.0, ti.f32)
        h = self.ps.support_radius
        # value of cubic spline smoothing kernel
        k = 1.0
        if self.ps.dim == 1:
            k = 4 / 3
        elif self.ps.dim == 2:
            k = 40 / 7 / np.pi
        elif self.ps.dim == 3:
            k = 8 / np.pi
        k /= h ** self.ps.dim
        q = r_norm / h
        if q <= 1.0:
            if q <= 0.5:
                q2 = q * q
                q3 = q2 * q
                res = k * (6.0 * q3 - 6.0 * q2 + 1)
            else:
                res = k * 2 * ti.pow(1 - q, 3.0)
        return res

    @ti.func
    def cubic_kernel_derivative(self, r):
        h = self.ps.support_radius
        # derivative of cubic spline smoothing kernel
        k = 1.0
        if self.ps.dim == 1:
            k = 4 / 3
        elif self.ps.dim == 2:
            k = 40 / 7 / np.pi
        elif self.ps.dim == 3:
            k = 8 / np.pi
        k = 6. * k / h ** self.ps.dim
        r_norm = r.norm()
        q = r_norm / h
        res = ti.Vector([0.0 for _ in range(self.ps.dim)])
        if r_norm > 1e-5 and q <= 1.0:
            grad_q = r / (r_norm * h)
            if q <= 0.5:
                res = k * q * (3.0 * q - 2.0) * grad_q
            else:
                factor = 1.0 - q
                res = k * (-factor * factor) * grad_q
        return res

    @ti.func
    def viscosity_force(self, p_i, p_j, r):
        # Compute the viscosity force contribution
        v_xy = (self.ps.v[p_i] -
                self.ps.v[p_j]).dot(r)
        res = 2 * (self.ps.dim + 2) * self.viscosity * (self.mass / (self.ps.density[p_j])) * v_xy / (
            r.norm()**2 + 0.01 * self.ps.support_radius**2) * self.cubic_kernel_derivative(
                r)
        return res

    @ti.func
    def pressure_force(self, p_i, p_j, r):
        # Compute the pressure force contribution, Symmetric Formula
        res = -self.density_0 * self.ps.m_V * (self.ps.pressure[p_i] / self.ps.density[p_i] ** 2
              + self.ps.pressure[p_j] / self.ps.density[p_j] ** 2) \
              * self.cubic_kernel_derivative(r)
        return res

    def substep(self):
        pass

    @ti.func
    def get_new_velocity(self, vn, vt) :
        mu_n = 0.98
        mu_t = 0.02

        a    = 1 - mu_t * (1 + mu_n) * vn.norm() / vt.norm()
        a    = ti.max(0,a)

        vn   = -mu_n * vn
        vt   = a     * vt

        return vn + vt

    @ti.func
    def segment_collision(self,p_i,q1,q2) :
        is_intersect = True
        dist_intersect = 0

        p1 = self.ps.x0[p_i]
        p2 = self.ps.x [p_i]

        p12 = p2-p1
        if p12.norm() <= 0 :
            is_intersect = False

        p12 = p12.normalized()
        q12 = (q2-q1).normalized()

        norm_p = ti.Vector([p12.y,-p12.x])
        norm_q = ti.Vector([q12.y,-q12.x])


        pa = norm_q.dot(p1-q1)
        pb = norm_q.dot(p2-q1)

        if pa * pb >= 0:
            if ti.abs(pb) >= self.ps.particle_radius :
                is_intersect = False

        qa = norm_p.dot(q1-p1)
        qb = norm_p.dot(q2-p1)

        if qa * qb >= 0:
            is_intersect = False

        if is_intersect :
            if pa < 0 :
                norm_q = - norm_q

            q1 += norm_q * self.ps.particle_radius 
            q1p2 = p2 - q1
            q1p2_on_q12 = q1p2.dot(q12) * q12
            
            self.ps.x[p_i] = q1 + q1p2_on_q12

            vec = (q1p2_on_q12 - q1p2).normalized()
            norm_vn_len = self.ps.v[p_i].dot(vec)
            if norm_vn_len < 0 :
                vn = norm_vn_len * vec
                vt = self.ps.v[p_i] - vn

                self.ps.v[p_i] = self.get_new_velocity(vn,vt)
        return is_intersect

    @ti.kernel
    def enforce_boundary(self):
        for p_i in range(self.ps.particle_num[None]):
            if self.ps.dim == 2:
                if self.ps.material[p_i] == self.ps.material_fluid:
                    for l_i in range(self.ps.lines.shape[0]//2) :
                        self.segment_collision(p_i,self.ps.lines[l_i*2+0],self.ps.lines[l_i*2+1]) 

    def step(self):
        self.ps.initialize_particle_system()
        self.substep()
        self.enforce_boundary()
