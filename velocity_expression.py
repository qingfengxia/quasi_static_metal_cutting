import math
import numpy as np
from parameter import *

if using_salome:
    from parameter_salome import *
else: 
    from parameter_gmsh import *

if workpiece_type_id  == 1:
    disc_H = 0.01;  #same with cutter now
    length_scale = disc_R;
    #if is_straight_chip:
    #    mesh_file = meshfolder + "/metal_cut_straight_chip.h5"
else:
    #rect work coord, already provided in salome
    length_scale = disc_R

if is_straight_chip:
    if using_salome:
        chip_Y = chip_top_y  # chip_length * math.cos(cutter_angle_v*math.pi/180);
    else:
        chip_Y = length_scale;
# still needed velocity expression

chip_speed = cutting_speed * (feed_thickness/chip_thickness)

frequency = cutting_speed / disc_R / 2 / math.pi
direction = -1
omega = cutting_speed / (disc_R+0.5*feed_thickness)

chip_omega = cutting_speed * (feed_thickness/chip_thickness)  / (chip_radius )  # correct
chip_sliding_vel_x = chip_speed*math.sin(cutter_angle_v*math.pi/180)
chip_sliding_vel_y = chip_speed*math.cos(cutter_angle_v*math.pi/180)
#print("chip_speed == chip_omega*(chip_radius+0.5*chip_thickness)?", chip_speed, chip_omega*(chip_radius+0.5*chip_thickness))


#chip_shear_angle = math.atan(-chip_start_y/chip_start_x);  # phi, recently renamed in gmsh

chip_friction_end_x = chip_friction_distance*math.sin(cutter_angle_v*math.pi/180);
chip_friction_end_y = chip_friction_distance*math.cos(cutter_angle_v*math.pi/180);
#chip_sliding_end_x = chip_sliding_distance*math.sin(cutter_angle_v*math.pi/180);
#chip_sliding_end_y = chip_sliding_distance*math.cos(cutter_angle_v*math.pi/180);

chip_center_x = chip_friction_end_x - (chip_radius+chip_thickness)*math.cos(cutter_angle_v*math.pi/180);
chip_center_y = chip_friction_end_y + (chip_radius+chip_thickness)*math.sin(cutter_angle_v*math.pi/180);

# need only by python code
p_chip_end_center_x = chip_center_x + (chip_radius + 0.5*chip_thickness) * math.cos(chip_end_angle*math.pi/180);
p_chip_end_center_y = chip_center_y + (chip_radius + 0.5*chip_thickness) * math.sin(chip_end_angle*math.pi/180);
p_chip_end_center_z = 0.5*cutter_thickness;
p_chip_end_xyz = p_chip_end_center_x, p_chip_end_center_y, p_chip_end_center_z
##################################################

if is_straight_chip:
    p_chip_end_o_x = chip_Y * math.tan(cutter_angle_v*math.pi/180);
    p_chip_end_o_y = chip_Y;
    p_chip_end_i_x = chip_Y * math.tan(cutter_angle_v*math.pi/180) - chip_thickness/math.cos(cutter_angle_v*math.pi/180);
    p_chip_end_i_y = chip_Y;
    def check_distance():
        #line p,q
        p = np.array([chip_friction_end_x, chip_friction_end_y, 0])
        q = np.array([p_chip_end_o_x, p_chip_end_o_y, 0])
        r1 = np.array([p_chip_end_i_x, p_chip_end_i_y, 0])
        r2 = np.array([chip_start_x, chip_start_y, 0])

        def t(p, q, r):
            x = p-q
            return np.dot(r-q, x)/np.dot(x, x)

        def d(p, q, r):
            return np.linalg.norm(t(p, q, r)*(p-q)+q-r)
        print('check distannce must match, ', d(p, q, r1), d(p, q, r2))
        # Prints 1.0
    check_distance()

if using_3D:
    dim = 3
else:
    dim = 2

###########################################
# Code for C++ evaluation of velocity
velocity_code = '''

class Velocity : public Expression
{
public:

  // Create expression with any components
  Velocity() : Expression(%d) {}

  // Function for evaluating expression on each cell
  void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
  {
    const double x0 = %f;
    const double y0 = %f;
    const double feed_thickness = %f;
    const double omega = %f;

    const uint cell_index = cell.index;
    const size_t did = (*subdomain_id)[cell_index];

    const double cx0 = %f;
    const double cy0 = %f;
    const double comega = %f;
    const double chip_sliding_end_x = %f;
    const double chip_sliding_end_y = %f;
    const double chip_sliding_vel_x = %f;
    const double chip_sliding_vel_y = %f;
    const double cutter_angle_v = %f;
    const double work_vel_x = y0*omega;

    const double turning_center_x = %f;
    const double turning_center_y = %f;
    const double fillet_r = %f;
    const double a_phi = %f;
    const double a_start = -pi/2;
    double a_delta = pi/2.0 - cutter_angle_v*pi/180;

    const int workpiece_type_id = %d;
    const int is_straight_chip = %d;
    const int using_fillet_shear_zone = %d;
    const int using_double_shear_heat_layer = %d;

    const int  work_subdomain_id = %d;
    const int  chip_subdomain_id = %d;
    const int  shear_subdomain_id = %d;
    const int  friction_subdomain_id = %d;

    values[0] = 0.0;
    values[1] = 0.0;
    //values[2] = 0.0;

    if(did == work_subdomain_id) {  // workpiece, left and right has diff radius an center
        if(workpiece_type_id == 1) {  // is disc
            double r = sqrt((x[0]-x0)*(x[0]-x0) + (x[1]-y0)*(x[1]-y0));
            double v = omega * r;
            double a = atan2((x[1]-y0), (x[0]-x0));
            if (x[0]<0) {
                double y0_1 = y0 + feed_thickness/2.0;
                r = sqrt((x[0]-x0)*(x[0]-x0) + (x[1]-y0_1)*(x[1]-y0_1));
                v = omega * r + feed_thickness/2.0;
                a = atan2((x[1]-y0_1), (x[0]-x0));
            }

            values[0] = -v * sin(a);
            values[1] = v * cos(a);
        }
        else { // workpiece rectangle
            values[0] = work_vel_x;  // only x-axis speed
            values[1] = 0.0;
        }
    }

    else if(did == chip_subdomain_id) {  // chip,  consisting of straight and arc sections
        if (is_straight_chip == 0) {
            double a = atan2((x[1]-cy0), (x[0]-cx0));
            //if (x[0] < chip_sliding_end_x && x[1] > chip_sliding_end_y) {
            if (a > (-cutter_angle_v*pi/180.0)) {
                double r = sqrt((x[0]-cx0)*(x[0]-cx0) + (x[1]-cy0)*(x[1]-cy0));
                double v = comega * r;
                values[0] = v * sin(a);
                values[1] = -v * cos(a);
            }
        }
        else {
            values[0] = chip_sliding_vel_x;
            values[1] = chip_sliding_vel_y;
        }
    }

    else if(did == shear_subdomain_id) {
        if(using_fillet_shear_zone) {// shear zone has the fillet 
            //double a_intersection = a_t - (pi - a_phi);
            double dist = sqrt((x[0]-turning_center_x)*(x[0]-turning_center_x) + (x[1]-turning_center_y)*(x[1]-turning_center_y));
            double shift = dist - fillet_r;
            double shifted_tcx = turning_center_x + shift * cos(a_phi);
            double shifted_tcy = turning_center_y - shift * sin(a_phi);
            double shifted_a = atan2((x[1]-shifted_tcy), (x[0]-shifted_tcx));
            double v_chip = sqrt(chip_sliding_vel_x*chip_sliding_vel_x + chip_sliding_vel_y*chip_sliding_vel_y);
            double v_feed = y0*omega;
            double v = v_chip + (1.0 - (shifted_a - a_start)/a_delta) * (v_feed - v_chip);
            values[0] = v * sin(-shifted_a);
            values[1] = v * cos(-shifted_a);
        }
        else {
            values[0] = chip_sliding_vel_x;
            values[1] = chip_sliding_vel_y;
            if(using_double_shear_heat_layer) {   // double mapping bc
                double a_t = atan2((x[1]), (x[0]));
                if ((a_t - (pi - a_phi)) > 1e-3) {
                    values[0] = work_vel_x;  // only x-axis speed
                    values[1] = 0;
                }
            }
        }
    }

    else if(did == friction_subdomain_id) {  // friction thin layer inside chip
        values[0] = chip_sliding_vel_x;
        values[1] = chip_sliding_vel_y;
    }

    else {
        values[0] = 0.0;
        values[1] = 0.0;
    }
  }
  // The data stored in mesh functions
  std::shared_ptr<MeshFunction<std::size_t> > subdomain_id;

};
'''%(dim, 0, -disc_R, feed_thickness, omega*direction, chip_center_x, chip_center_y, chip_omega*direction,
        chip_friction_end_x, chip_friction_end_y, chip_sliding_vel_x, chip_sliding_vel_y, cutter_angle_v,
        turning_center_x, turning_center_y, fillet_r,  shear_angle *  math.pi/180, 
        workpiece_type_id, is_straight_chip, int(using_fillet_shear_zone), int(using_double_shear_heat_layer),
        work_subdomain_id, chip_subdomain_id, shear_subdomain_id, friction_subdomain_id
        )

#print(velocity_code)