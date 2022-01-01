import taichi as ti
import numpy as np
from particle_system import ParticleSystem
from wcsph import WCSPHSolver

ti.init(arch=ti.cpu)
# ti.init(arch=ti.cpu,kernel_profiler=True)

# Use GPU for higher peformance if available
# ti.init(arch=ti.gpu, device_memory_GB=1, packed=True, default_fp=ti.f32)
# ti.init(arch=ti.gpu, device_memory_GB=1.1, packed=True, kernel_profiler=True)

if __name__ == "__main__":
    ps = ParticleSystem((512, 512))
    # ps = ParticleSystem((256, 256))

    ps.add_cube(lower_corner=[5, 1],
                cube_size=[4.0, 3.0],
                velocity=[0.0, 0.0],
                density=1000.0,
                color=0x956333,
                material=1)

    ps.add_cube(lower_corner=[1, 1],
                cube_size=[4.0, 3.0],
                velocity=[0.0, 0.0],
                density=1000.0,
                color=0x956333,
                material=1)

    wcsph_solver = WCSPHSolver(ps)
    gui = ti.GUI('WCSPH',ps.res,background_color=0xFFFFFF)
    gui.fps_limit = 20

    is_mov = 0
    while gui.running:

        for e in gui.get_events(ti.GUI.PRESS):
            if e.key == ti.GUI.LMB:
                is_mov = 1 - is_mov

        if is_mov == 1 :
            epos = np.array(gui.get_cursor_pos())
            epos = epos / ps.screen_to_world_ratio * ps.res[0]
            ps.mv_cross(epos[0], epos[1])

        for i in range(5):
            wcsph_solver.step()

        particle_info = ps.dump()
        num   = len(particle_info['position'])

        part1 = particle_info['position'][:num//2,:]
        part2 = particle_info['position'][num//2:,:]

        gui.circles(part1 * ps.screen_to_world_ratio / ps.res[0],
                    radius=ps.particle_radius / 1.5 * ps.screen_to_world_ratio,
                    color=0x0066ff)
        gui.circles(part2 * ps.screen_to_world_ratio / ps.res[0],
                    radius=ps.particle_radius / 1.5 * ps.screen_to_world_ratio,
                    color=0x33ffff)

        for i in range(particle_info['lines'].shape[0]//2) :
            p1 = particle_info['lines'][i*2+0,:] * ps.screen_to_world_ratio / ps.res[0]
            p2 = particle_info['lines'][i*2+1,:] * ps.screen_to_world_ratio / ps.res[0]
            gui.line(begin=p1, end=p2, radius=2, color=0x7f7f7f)

        gui.show()
    
    # ti.print_kernel_profile_info()
