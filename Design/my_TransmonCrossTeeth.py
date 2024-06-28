### Haimeng

from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import BaseQubit
import numpy as np
class myTransmonNewCrosswTeeth(BaseQubit):
    
    """Component metadata"""
    default_options = Dict(
        make_teeth=True,
        cross_type_reduced = False,
        pin_orientation = 0,
        flux_orientation = 90,
        north_width='10um',
        north_gap='6um',
        cap_width='10um',
        cap_gap='6um',
        cap_gap_ground='6um',
        finger_length='20um',
        finger_count='5',
        cap_distance='50um',
        
       
        south_width='10um',
        south_gap='6um',
        termination_cap_width='10um',
        termination_cap_gap='6um',
        termination_cap_gap_ground='6um',
        termination_finger_length='20um',
        termination_finger_count='5',
        termination_cap_distance='50um',
        
        cross_width='20um',
        cross_length='200um',
        cross_gap='20um',
        chip='main',
        _default_connection_pads=Dict(
        connector_type='0',  # 0 = Claw type, 1 = gap type
        claw_length='30um',
        ground_spacing='5um',
        claw_width='10um',
        claw_gap='6um',
        connector_location=
        '0'  # 0 => 'west' arm, 90 => 'north' arm, 180 => 'east' arm
        ),
        make_fl = True,
        fl_style=None,
        fl_options=Dict(t_top='15um',
                        t_offset='0um',
                        t_inductive_gap='3um',
                        t_width='5um',
                        t_gap='3um'),
        taper_options=Dict(t_width_i='10um',
                            t_width_f='4.25um', 
                            t_length='90um', 
                            t_gap_i ='6um', 
                            t_gap_f ='2.5um', 
                            t_punch_through = '5um',
                            t_hanger_para_length='40um',
                            t_hanger_para_width='2.5um',
                            t_hanger_para_gap = '2.5um',
                            t_hanger_perp_length='55um',
                            t_hanger_perp_width='2.7um',
                            t_hanger_perp_gap = '3.2um')

    )
        
    component_metadata = Dict(short_name='Cross',
                          _qgeometry_table_poly='True',
                          _qgeometry_table_junction='True')
    
    def make(self):
        self.make_pocket()
        self.make_connection_pads()
        if self.options.make_fl == True:
            self.make_flux_line_tapered()
    
    def make_pocket(self):
        
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p

        cross_width = p.cross_width
        cross_length = p.cross_length
        cross_gap = p.cross_gap

        # access to chip name
        chip = p.chip
#         p.finger_length = int(p.finger_count)/2 * (float(p.cap_width) + float(p.cap_gap)) + 0.5 * float(p.north_width) -float(p.cap_width) 
        # Creates the cross and the etch equivalent.
        cross_type_reduced = p.cross_type_reduced
        if cross_type_reduced:
            cross_line = draw.shapely.ops.unary_union([
                draw.LineString([(0, cross_length), (0, -cross_length)]),
                draw.LineString([(cross_length, 0), (0, 0)])
            ])
        else:
                        cross_line = draw.shapely.ops.unary_union([
                draw.LineString([(0, cross_length), (0, -cross_length)]),
                draw.LineString([(cross_length, 0), (-cross_length, 0)])
            ])

        cross = cross_line.buffer(cross_width / 2, cap_style=2)
        cross_etch = cross.buffer(cross_gap, cap_style=3, join_style=2)

        # The junction/SQUID
        #rect_jj = draw.rectangle(cross_width, cross_gap)
        #rect_jj = draw.translate(rect_jj, 0, -cross_length-cross_gap/2)
        rect_jj = draw.LineString([(cross_length,0),
                                   (cross_length + cross_gap, 0)])

      
        
        #### NORTH TEETH
        name = 'teeth'
        p = self.p
        N = int(p.finger_count)

        #Finger Capacitor
        cap_box = draw.rectangle(N * p.cap_width + (N - 1) * p.cap_gap,
                                 p.cap_gap + 2 * p.cap_width + p.finger_length,
                                 0, 0)
        
        make_cut_list = []
        make_cut_list.append([0, (p.finger_length) / 2])
        make_cut_list.append([(p.cap_width) + (p.cap_gap / 2),
                              (p.finger_length) / 2])
        flip = -1

        for i in range(1, N):
            make_cut_list.append([
                i * (p.cap_width) + (2 * i - 1) * (p.cap_gap / 2),
                flip * (p.finger_length) / 2
            ])
            make_cut_list.append([
                (i + 1) * (p.cap_width) + (2 * i + 1) * (p.cap_gap / 2),
                flip * (p.finger_length) / 2
            ])
            flip = flip * -1

        cap_cut = draw.LineString(make_cut_list).buffer(p.cap_gap / 2,
                                                        cap_style=2,
                                                        join_style=2)
        cap_cut = draw.translate(cap_cut,
                                 -(N * p.cap_width + (N - 1) * p.cap_gap) / 2,
                                 0)

        cap_body = draw.subtract(cap_box, cap_cut)
        cap_body = draw.translate(
            cap_body, 0, -p.cap_distance -
            (p.cap_gap + 2 * p.cap_width + p.finger_length) / 2)

        cap_etch = draw.rectangle(
            N * p.cap_width + (N - 1) * p.cap_gap + 2 * p.cap_gap_ground,
            p.cap_gap + 2 * p.cap_width + p.finger_length +
            2 * p.cap_gap_ground, 0, -p.cap_distance -
            (p.cap_gap + 2 * p.cap_width + p.finger_length) / 2)

        #CPW
        north_cpw = draw.LineString([[0, 0], [0, -p.cap_distance]])
        
         #Rotate and Translate
        c_items = [north_cpw, cap_body, cap_etch]
        
        c_items = draw.translate(c_items, 0,  2* p.cap_width + p.cap_distance + p.cap_gap + p.finger_length + p.cross_length )
        [north_cpw, cap_body, cap_etch] = c_items
      
        
        north_pin_list = north_cpw.coords
        

        self.add_pin('res_pin',
                     points=np.array(north_pin_list[::-1]),
                     width=p.north_width,
                     input_as_norm=True)
        
        
        ### SOUTH TEETH (TERMINATION TEETH)
        
        name = 'termination'
        p = self.p
        N = int(p.termination_finger_count)

        #Finger Capacitor
        termination_cap_box = draw.rectangle(N * p.termination_cap_width + (N - 1) * p.termination_cap_gap,
                                 p.termination_cap_gap + 2 * p.termination_cap_width + p.termination_finger_length,
                                 0, 0)
        
        termination_make_cut_list = []
        termination_make_cut_list.append([0, (p.termination_finger_length) / 2])
        termination_make_cut_list.append([(p.termination_cap_width) + (p.termination_cap_gap / 2),
                              (p.termination_finger_length) / 2])
        flip = -1

        for i in range(1, N):
            termination_make_cut_list.append([
                i * (p.termination_cap_width) + (2 * i - 1) * (p.termination_cap_gap / 2),
                flip * (p.termination_finger_length) / 2
            ])
            termination_make_cut_list.append([
                (i + 1) * (p.termination_cap_width) + (2 * i + 1) * (p.termination_cap_gap / 2),
                flip * (p.termination_finger_length) / 2
            ])
            flip = flip * -1

        termination_cap_cut = draw.LineString(termination_make_cut_list).buffer(p.termination_cap_gap / 2,
                                                        cap_style=2,
                                                        join_style=2)
        termination_cap_cut = draw.translate(termination_cap_cut,
                                 -(N * p.termination_cap_width + (N - 1) * p.termination_cap_gap) / 2,
                                 0)

        termination_cap_body = draw.subtract(termination_cap_box, termination_cap_cut)
        termination_cap_body = draw.translate(
            termination_cap_body, 0, -(-p.termination_cap_distance -
            (p.termination_cap_gap + 2 * p.termination_cap_width + p.termination_finger_length) / 2))

        termination_cap_etch = draw.rectangle(
            N * p.termination_cap_width + (N - 1) * p.termination_cap_gap + 2 * p.termination_cap_gap_ground,
            
            p.termination_cap_gap + 2 * p.termination_cap_width + p.termination_finger_length +
            2 * p.termination_cap_gap_ground, 
            0, 
            -(-p.termination_cap_distance -
            (p.termination_cap_gap + 2 * p.termination_cap_width + p.termination_finger_length) / 2))

        #CPW
        south_cpw = draw.LineString([[0, 0], [0, p.termination_cap_distance]])
        
        #Rotate and Translate
        termination_c_items = [south_cpw, termination_cap_body, termination_cap_etch, termination_cap_cut]
        termination_c_items = draw.translate(termination_c_items, 0, -(  2* p.termination_cap_width + p.termination_cap_distance + p.termination_cap_gap + p.termination_finger_length + p.cross_length ))
        [south_cpw, termination_cap_body, termination_cap_etch, termination_cap_cut] = termination_c_items
        
        # Cross
        #rotate and translate
        polys = [cross, cross_etch, rect_jj]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)

        [cross, cross_etch, rect_jj] = polys
        
        pins = [cap_body, termination_cap_body, north_cpw, south_cpw, cap_etch,termination_cap_etch ]
        pins = draw.rotate(pins, p.pin_orientation, origin=(0, 0))
        pins = draw.translate(pins, p.pos_x, p.pos_y)

        [cap_body, termination_cap_body, north_cpw, south_cpw, cap_etch,termination_cap_etch ] = pins
        
        
        
          # generate qgeometry
        cross_combined = draw.union(cross, cap_body, termination_cap_body)
        self.add_qgeometry('poly', dict(cross=cross_combined), chip=chip)
        self.add_qgeometry('poly',
                           dict(cross_etch=cross_etch),
                           subtract=True,
                           chip=chip)
        self.add_qgeometry('junction',
                           dict(rect_jj=rect_jj),
                           width=cross_width,
                           chip=chip)
        
        self.add_qgeometry('path', {'north_cpw': north_cpw},
                           width=p.north_width,
                           layer=p.layer)
        self.add_qgeometry('path', {'north_cpw_sub': north_cpw},
                           width=p.north_width + 2 * p.north_gap,
                           layer=p.layer,
                           subtract=True)
        
                
        self.add_qgeometry('path', {'south_cpw': south_cpw},
                           width=p.south_width,
                           layer=p.layer)
        self.add_qgeometry('path', {'south_cpw_sub': south_cpw},
                           width=p.south_width + 2 * p.south_gap,
                           layer=p.layer,
                           subtract=True)
        

        #self.add_qgeometry('poly', {'cap_body': cap_body}, layer=p.layer)
        self.add_qgeometry('poly', {'cap_etch': cap_etch},
                           layer=p.layer,
                           subtract=True)
        self.add_qgeometry('poly', {'termination_cap_etch': termination_cap_etch},
                   layer=p.layer,
                   subtract=True)
        

#         points = list(draw.shapely.geometry.shape(north_cpw).coords)
#         self.add_pin(name, points, p.north_width)  # TODO: chip
        
        north_pin_list = north_cpw.coords
        

        self.add_pin('res_pin',
                     points=np.array(north_pin_list[::-1]),
                     width=p.north_width,
                     input_as_norm=True)
        
        south_pin_list = south_cpw.coords
        

        self.add_pin('ground_pin',
                     points=np.array(south_pin_list[::-1]),
                     width=p.south_width,
                     input_as_norm=True)
        
        
    def make_connection_pads(self):
        """Goes through connector pads and makes each one."""
        for name in self.options.connection_pads:
            self.make_connection_pad(name)


    def make_connection_pad(self, name: str):
        """Makes individual connector pad.

        Args:
            name (str) : Name of the connector pad
        """

        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        cross_width = p.cross_width
        cross_length = p.cross_length
        cross_gap = p.cross_gap
        
        # access to chip name
        chip = p.chip

        pc = self.p.connection_pads[name]  # parser on connector options
        c_g = pc.claw_gap
        c_l = pc.claw_length
        c_w = pc.claw_width
        g_s = pc.ground_spacing
        con_loc = pc.connector_location

        claw_cpw = draw.box(0, -c_w / 2, -4 * c_w, c_w / 2)

        if pc.connector_type == 0:  # Claw connector
            t_claw_height = 2*c_g + 2 * c_w + 2*g_s + \
                2*cross_gap + cross_width  # temp value

            claw_base = draw.box(-c_w, -(t_claw_height) / 2, c_l,
                                 t_claw_height / 2)
            claw_subtract = draw.box(0, -t_claw_height / 2 + c_w, c_l,
                                     t_claw_height / 2 - c_w)
            claw_base = claw_base.difference(claw_subtract)

            connector_arm = draw.shapely.ops.unary_union([claw_base, claw_cpw])
            connector_etcher = draw.buffer(connector_arm, c_g)
        else:
            connector_arm = claw_cpw
            connector_etcher = draw.buffer(connector_arm, c_g)

        # Making the pin for  tracking (for easy connect functions).
        # Done here so as to have the same translations and rotations as the connector. Could
        # extract from the connector later, but since allowing different connector types,
        # this seems more straightforward.
        # port_line = draw.LineString([(-4 * c_w, -c_w / 2), (-4 * c_w, c_w / 2)])

        claw_rotate = 0
        if con_loc > 135:
            claw_rotate = 180
        elif con_loc > 45:
            claw_rotate = -90

        # Rotates and translates the connector polygons (and temporary port_line)
        polys = [connector_arm, connector_etcher, port_line]
        polys = draw.translate(polys, -(cross_length + cross_gap + g_s + c_g),
                               0)
        polys = draw.rotate(polys, claw_rotate, origin=(0, 0))
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [connector_arm, connector_etcher, port_line] = polys

        # Generates qgeometry for the connector pads
        self.add_qgeometry('poly', {f'{name}_connector_arm': connector_arm},
                           chip=chip)
        self.add_qgeometry('poly',
                           {f'{name}_connector_etcher': connector_etcher},
                           subtract=True,
                           chip=chip)

        # self.add_pin(name, port_line.coords, c_w)

    def make_flux_line(self):
        """Creates the charge line if the user has charge line option to
        TRUE."""

        # Grab option values
        pf = self.p.fl_options
        p = self.p
        #Make the T flux line
        h_line = draw.LineString([(-pf.t_top / 2, 0), (pf.t_top / 2, 0)])
        v_line = draw.LineString([(pf.t_offset, 0), (pf.t_offset, -0.03)])

        parts = [h_line, v_line]

        # Move the flux line down to the SQUID
        parts = draw.translate(
            parts, 0, -(p.cross_length + p.cross_gap + pf.t_inductive_gap +
                        pf.t_width / 2 + pf.t_gap))

        # Rotate and translate based on crossmon location
        parts = draw.rotate(parts, p.flux_orientation, origin=(0, 0))
        parts = draw.translate(parts, p.pos_x, p.pos_y)

        [h_line, v_line] = parts

        # Adding to qgeometry table
        self.add_qgeometry('path', {
            'h_line': h_line,
            'v_line': v_line
        },
                           width=pf.t_width,
                           layer=p.layer)

        self.add_qgeometry('path', {
            'h_line_sub': h_line,
            'v_line_sub': v_line
        },
                           width=pf.t_width + 2 * pf.t_gap,
                           subtract=True,
                           layer=p.layer)

        # Generating pin
        pin_line = v_line.coords
        self.add_pin("flux_line",
                     points=pin_line,
                     width=pf.t_width,
                     gap=pf.t_gap,
                     input_as_norm=True)

    def make_flux_line_tapered(self):
        """
        Creates the charge line which is taper w/ a hanging bar.
        
        Activates when all of these conditions are met:
        - self.options.make_fl == True
        - self.options.fl_options.ll_style == "tapered"
        """
        # Grab option values
        pf = self.p.taper_options
        p = self.p
        
        
        ##### Fast Flux Line Element #####
        v_line = draw.LineString([(-0.0, pf.t_width_f/4), (0.0, -pf.t_width_f/4)])
        # sub_v_line = draw.LineString([(0.0, pf.t_width_f/2 + pf.t_hanger_perp_gap), (0.01, - pf.t_width_f/2 - pf.t_hanger_perp_gap)])
        # Make tappered FF Line
        tapper = draw.Polygon([(pf.t_width_i / 2, -pf.t_width_i / 4), 
                               (- pf.t_width_i/2, -pf.t_width_i / 4), 
                               (-pf.t_width_f / 2, pf.t_length), 
                               (-pf.t_width_f/2 + pf.t_hanger_para_length, pf.t_length),
                               (-pf.t_width_f/2 + pf.t_hanger_para_length, pf.t_length - pf.t_hanger_perp_length),
                               (-pf.t_width_f/2 + pf.t_hanger_para_length - pf.t_hanger_perp_width, pf.t_length - pf.t_hanger_perp_length),
                               (-pf.t_width_f/2 + pf.t_hanger_para_length - pf.t_hanger_perp_width, pf.t_length - pf.t_hanger_para_width),
                               ((pf.t_length - pf.t_hanger_para_width)*(pf.t_width_f - pf.t_width_i) / (2 * pf.t_length) + pf.t_width_i/2, pf.t_length - pf.t_hanger_para_width)])
        
        # Make subtraction      
        sub_tapper = draw.Polygon([(pf.t_width_i / 2 + pf.t_gap_i, -pf.t_width_i / 4), 
                                   (- pf.t_width_i/2 - pf.t_gap_i, -pf.t_width_i / 4), 
                                   (-pf.t_width_f / 2 - pf.t_gap_f, pf.t_length + pf.t_hanger_para_gap), 
                                   (-pf.t_width_f/2 + pf.t_hanger_para_length + pf.t_hanger_perp_gap, pf.t_length + pf.t_hanger_para_gap),
                                   (-pf.t_width_f/2 + pf.t_hanger_para_length + pf.t_hanger_perp_gap, pf.t_length - pf.t_hanger_perp_length),
                                   (-pf.t_width_f/2 + pf.t_hanger_para_length - pf.t_hanger_perp_width - pf.t_hanger_perp_gap, pf.t_length - pf.t_hanger_perp_length),
                                   (-pf.t_width_f/2 + pf.t_hanger_para_length - pf.t_hanger_perp_width - pf.t_hanger_perp_gap, pf.t_length - pf.t_hanger_para_width - pf.t_hanger_para_gap),
                                   ((pf.t_length - pf.t_hanger_para_width - pf.t_hanger_para_gap)/((pf.t_length - pf.t_hanger_para_gap)/((pf.t_width_f - pf.t_width_i)/2 + pf.t_gap_f - pf.t_gap_i)) + pf.t_width_i / 2 + pf.t_gap_i,pf.t_length - pf.t_hanger_para_width - pf.t_hanger_para_gap)])
        
        # Translate all the parts to the edge of the 
        parts = [v_line,tapper, sub_tapper]
    
        parts = draw.translate(
            parts, 0, -(p.cross_length + p.cross_gap + pf.t_length - pf.t_punch_through))
        parts = draw.rotate(parts, p.flux_orientation, origin=(0, 0))
        parts = draw.translate(parts, p.pos_x, p.pos_y)
        
        v_line,tapper, sub_tapper = parts
        
        
        
        ##### Adding to qgeometry table #####
        # self.add_qgeometry('poly', dict(v_line=v_line), layer=p.layer)
        # self.add_qgeometry('poly', dict(v_line=sub_v_line), subtract=True, layer=p.layer)
        self.add_qgeometry('poly', dict(tapper=tapper), layer=p.layer)
        self.add_qgeometry('poly', 
                           dict(sub_tapper=sub_tapper), 
                           subtract=True, 
                           layer=p.layer)
         # Generating pin
        pin_line = v_line.coords
        self.add_pin("flux_line",
                     points=pin_line,
                     width=pf.t_width_f,
                     gap=pf.t_gap,
                     input_as_norm=True)
            