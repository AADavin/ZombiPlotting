import ete3 
import drawSvg as draw
import random

class SpeciesTreeDrawer():
    def __init__(self, treefile, eventsfile, height, width, padding=15, background=None):          
        
        self.height = height        
        self.width =  width
        self.padding = padding        
        self.branch2features, self.n_leaves, self.total_time = self._read_eventsfile(eventsfile)
        self.tree = self._read_treefile(treefile)                
        self.node2coor = self._get_coordinates()                        
        self.d = draw.Drawing(self.width, self.height, origin=(0,0), displayInline=False) 
        
        if background != None:
            self.d.draw(draw.Rectangle(0,0,width,height, fill=background))

    def _read_eventsfile(self, eventsfile):
        
        n_leaves = 0
        branch2features = dict()
        branch2features["Root"] = {
                        "ori_time": 0,
                        "end_time": 0,
                        "type": "O"}        
        with open(eventsfile) as f:
            f.readline()
            for l in f:
                time, event, nodes = l.strip().split("\t")                       
                if event == "S":
                    node_p, node_c1, node_c2 = nodes.split(";")
                    branch2features[node_p]["end_time"] = float(time)
                    branch2features[node_c1] = dict()
                    branch2features[node_c2] = dict()
                    branch2features[node_c1]["ori_time"] = float(time)
                    branch2features[node_c2]["ori_time"] = float(time)            
                elif event == "E":
                    branch2features[nodes]["end_time"] = float(time)
                    n_leaves +=1
                elif event == "F":
                    branch2features[nodes]["end_time"] = float(time)
                    n_leaves +=1
                    total_time = float(time)
        return branch2features, n_leaves, total_time
    
            
    def _read_treefile(self, treefile):
        with open(treefile) as f:    
            tree = ete3.Tree(f.readline().strip(), format=1)
        return(tree)
    
    def _y_scaling(self, x):    
        return round(self.height * x / self.total_time, 2)
    
    def _get_coordinates(self):    
        
        node2coor = dict()
        padding = self.padding
        spacing = (self.width - padding) / (self.n_leaves + 1)
        total_time = self.total_time
        leaf_counter = 0

        for n in self.tree.traverse("postorder"):    

            dist_to_origin = n.get_distance(self.tree) + self.tree.dist    

            if n.is_leaf():
                leaf_counter +=1
                x_coor = padding + (spacing * leaf_counter)
            else:        
                c1, c2 = n.get_children()
                c1_xcoor = c1.x_coor
                c2_xcoor = c2.x_coor
                mean_point = abs(c1_xcoor - c2_xcoor) / 2.0
                x_coor = min([c1_xcoor, c2_xcoor]) + mean_point
            
            y_coor = self._y_scaling(total_time) - self._y_scaling(dist_to_origin) 
            n.add_feature("x_coor", round(x_coor,2))
            n.add_feature("y_coor", round(y_coor,2))  
            node2coor[n.name] = (n.x_coor, n.y_coor)
        return node2coor
    
    def DrawCells(self, url, stroke="red", stroke_width=1, radius = 3, inter_distance = 2, fill="transparent", draw_h = True, h_line_width=1, draw_out = False):   
        
        spacing = self.width / (self.n_leaves + 1)
        width = self.width
        height = self.height
        d = self.d
           
        #d.append(draw.Circle(150, 100, radius, stroke_width=2, stroke='white', stroke_dasharray="5,5", fill="transparent"))
        #d.append(draw.Ellipse(150, 100, 80, 60, stroke_width=3, stroke='white', fill ="transparent"))
       
        for n in self.tree.traverse():
            # Inner tree, vertical lines
            x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)

            n_cells = int(y_1 / ((radius * 2) + inter_distance ))
            
            for i in range(n_cells):            
                d.append(draw.Circle(x_0, y_0 + (inter_distance * i) + (i * radius * 2) + inter_distance, radius, stroke_width=stroke_width, stroke=stroke, fill=fill))
                if draw_out:
                    if i == 0:
                        continue
                    p = draw.Path(stroke_width=stroke_width, stroke=stroke, fill="none")
                    p.M(x_0, y_0 + (inter_distance * i) + (i * radius * 2) + inter_distance).l(-radius*3, -radius*3)
                    d.append(p)
            
            if draw_h:
                if not n.is_leaf():        
                    c1, c2 = n.get_children()                
                    # Horizontal lines        
                    # Inner Tree                
                    p = draw.Path(stroke=stroke, stroke_width=h_line_width, fill='none')  
                    p.M(c1.x_coor, n.y_coor).L(c2.x_coor, n.y_coor)
                    d.append(p)      

            
        d.setRenderSize(h=self.height)
        d.saveSvg(url)
        return(d)
            
    def DrawInner(self, url, stroke="red", stroke_width=1):   
        
        spacing = self.width / (self.n_leaves + 1)
        total_time = self.total_time
        width = self.width
        height = self.height
        d = self.d
        stroke_a = stroke_width / 2.0
       
        for n in self.tree.traverse():
            # Inner tree, vertical lines
            x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)
            p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
            p.M(x_0, y_0).l(0, y_1)  
            d.append(p)                        
            if not n.is_leaf():        
                c1, c2 = n.get_children()                
                # Horizontal lines        
                # Inner Tree                
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(c1.x_coor - stroke_a, n.y_coor).L(c2.x_coor + stroke_a, n.y_coor)
                d.append(p)            
        d.setRenderSize(h=self.height)
        d.saveSvg(url)
        return(d)
    
    def add_names(self, url, font_size = 12, family = "Arial", colour="white", x_adj = 0, y_adj = 0):
        d = self.d
        for n in self.tree.traverse():
            x_0, y_0 = n.x_coor, n.y_coor
            p = draw.Text(n.name,font_size, x_0 + x_adj, y_0 + y_adj, center=False, fill=colour, font_family=family)
            d.append(p)
        d.saveSvg(url)
        return(d)
    
    def DrawOuter(self, url, stroke="white", stroke_width=1, radius = 3):
            
        spacing = self.width / (self.n_leaves + 1)        
        total_time = self.total_time
        width = self.width
        height = self.height
        d = self.d
        stroke_a = stroke_width / 2.0
    
        for n in self.tree.traverse():   
            
            if n.is_leaf():

                upx, upy = (n.up).x_coor, (n.up).y_coor 
                x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist - radius)
                
                if x_0 < upx:          
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) - radius
                    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                    p.M(x_0 + radius, y_0).l(0, y_1)
                    d.append(p) 
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) + radius
                    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                    p.M(x_0 - radius, y_0).l(0, y_1)
                    d.append(p)                     
                else:                    
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) + radius
                    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                    p.M(x_0 + radius, y_0).l(0, y_1)
                    d.append(p)                 
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) - radius
                    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                    p.M(x_0 - radius, y_0).l(0, y_1)
                    d.append(p) 
            else:    
                
                       
                x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)  

                # Horizontal lines

                c1, c2 = n.get_children()

                # Line below
                x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)                    
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(c1.x_coor + radius - stroke_a, n.y_coor - radius).L(c2.x_coor - radius + stroke_a, n.y_coor - radius)
                d.append(p)                       


                x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)                    
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(c1.x_coor - radius - stroke_a, n.y_coor + radius).L(n.x_coor - radius + stroke_a, n.y_coor + radius)
                d.append(p) 

                x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)                    
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(n.x_coor + radius - stroke_a, n.y_coor + radius).L(c2.x_coor + radius + stroke_a, n.y_coor + radius)
                d.append(p) 
                    
                    
                if not n.is_root():          
                                 
                    # Vertical lines
                    
                    upx, upy = (n.up).x_coor, (n.up).y_coor    
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) 
                    
                    
                    if x_0 < upx:          
                        x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) 
                        p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                        p.M(x_0 - radius, y_0 + radius).l(0, y_1)
                        d.append(p) 
                        x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) 
                        p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                        p.M(x_0 + radius, y_0 + radius).l(0, y_1 - (radius * 2))
                        d.append(p)                     
                    else:                    
                        x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) 
                        p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                        p.M(x_0 - radius, y_0 + radius).l(0, y_1 - (radius*2))
                        d.append(p)                 
                        x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) 
                        p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                        p.M(x_0 + radius, y_0 + radius).l(0, y_1)
                        d.append(p) 
                else:
                    
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) 
                    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                    p.M(x_0 - radius, y_0 + radius).l(0, y_1)
                    d.append(p) 
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist) 
                    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                    p.M(x_0 + radius, y_0 + radius).l(0, y_1)
                    d.append(p)                                             

        d.setRenderSize(h=self.height)
        d.saveSvg(url)
        return(d)
    
    
    def DrawGeneTree(self, url, genetreefile, eventsgenetreefile, stroke="red", stroke_width=2, displacement=2, r_displacement = 0, symbol_scaling=1.5):
        
        spacing = self.width / (self.n_leaves + 1)
        total_time = self.total_time
        width = self.width
        height = self.height
        d = self.d       
        radius = 0.1
        
        genetree_svg = list()
        
        r_displacement = random.uniform(-r_displacement, r_displacement)
        node2events = dict()
        node2displacement = dict()

        
        gtree = self._read_treefile(genetreefile)       
        
        
        #################### This a patch to correct a bug in Zombi ################
        
        with open(eventsgenetreefile) as f:
            f.readline()        
            o, e, nodes = f.readline().strip().split("\t")
            
            o_time = float(o)
            while e != "S" and e != "E" and e !="F" and e != "L" and e !="T" and e !="D":
                t, e, nodes = f.readline().strip().split("\t")
            gtree.dist = float(t) - float(o)            

          
        #############################################################################
        
        # We add the coordinates
        for n in gtree.traverse():            
            dist_to_origin = n.get_distance(gtree) + gtree.dist + o_time
            x_coor, _ = self.node2coor[n.name.split("_")[0]] 
            y_coor = round(self._y_scaling(total_time) - self._y_scaling(dist_to_origin),2)  
            
            n.add_feature("displacement", 0)
            n.add_feature("x_coor", x_coor)
            n.add_feature("y_coor", y_coor)    
            
        # We add the displacement
        gtree.displacement = r_displacement
        gtree.x_coor += r_displacement
        for n in gtree.iter_descendants():
            parent_d = (n.up).displacement
            n.displacement = parent_d
            if n.name in node2displacement:
                n.displacement += node2displacement[n.name]
            n.x_coor += (n.displacement*displacement) + r_displacement
            #n.y_coor += (n.displacement*4)
            
        # We get a dict with the new coordinates in the gene tree
        
        gnode2x_coor = {n.name:n.x_coor for n in gtree.traverse()}
        
        
        ###### EVENT SYMBOLS ############
       
        arrow = draw.Marker(-0.1, -0.5, 0.9, 0.5, scale=stroke_width*1.5, orient='auto')
        arrow.append(draw.Lines(-0.1, -0.5, -0.1, 0.5, 0.9, 0, fill=stroke, close=True))
 
        ##################################
    
         # We add the symbols
        symbols = list()
        
        with open(eventsgenetreefile) as f:
            f.readline()
            for l in f:
                t, e, nodes = l.strip().split("\t")                
                if nodes == "Root":
                    x_coor = gnode2x_coor["Root_1"]   
                else:
                    x_coor = gnode2x_coor["_".join(nodes.split(";")[:2])]   
                if e == "O":
                    o_time = float(t)                      
                    y_coor = self._y_scaling(total_time) - self._y_scaling(float(t))
                    symbols.append(("O", x_coor, y_coor)) 
                if e == "D":
                    _, _, node1, num1, node2, num2 = nodes.split(";")
                    n1 = node1 + "_" + num1
                    n2 = node2 + "_" + num2
                    node2displacement[n1] = 1
                    node2displacement[n2] = -1  
                    y_coor = self._y_scaling(total_time) - self._y_scaling(float(t))
                    symbols.append(("D", x_coor, y_coor))
                if e == "L":
                    y_coor = self._y_scaling(total_time) - self._y_scaling(float(t))
                    symbols.append(("L", x_coor, y_coor)) 
                if e == "T":   
                    _, _, node1, num1, node2, num2 = nodes.split(";")
                    n1 = node1 + "_" + num1
                    n2 = node2 + "_" + num2
                    node2events["_".join(nodes.split(";")[0:2])] = "T"            
                    node2displacement[n2] = 1
                    #x_coor, _ = self.node2coor[nodes.split(";")[4]]
                    #y_coor = self._y_scaling(float(t))   
 
        for n in gtree.traverse():
            
            # Need to check 2 things: is leaf? and has_events?
            
            if n.is_leaf() and n.name in node2events:
                
                # Vertical line
                
                if node2events[n.name] == "L":
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)
                    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none', marker_end=cross)  
                    p.M(x_0, y_0).l(0, y_1)  
                else:
                    x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)
                    p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                    p.M(x_0, y_0).l(0, y_1)  
                genetree_svg.append(p)

                
            elif not n.is_leaf() and n.name in node2events:
                
                # Vertical line
                
                x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(x_0, y_0).l(0, y_1)
                genetree_svg.append(p)
                
                # Horizontal line
                
                c1, c2 = n.get_children() 
                
                if node2events[n.name] == "T":                    
                        p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none', 
                                      stroke_dasharray="4", marker_end=arrow)  

                        if c1.x_coor > c2.x_coor:                    
                            p.M(c1.x_coor, n.y_coor).L(c2.x_coor + radius, n.y_coor)
                        else:
                            p.M(c1.x_coor, n.y_coor).L(c2.x_coor - radius, n.y_coor)

                genetree_svg.append(p)
                
            
            elif n.is_leaf() and n.name not in node2events:
                
                # Vertical line
                
                x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(x_0, y_0).l(0, y_1)  
                genetree_svg.append(p)

                
            elif not n.is_leaf() and n.name not in node2events:
                
                # Vertical line
                
                x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(x_0, y_0).l(0, y_1)  
                genetree_svg.append(p)

                
                # Horizontal line
                
                c1, c2 = n.get_children() 
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(c1.x_coor - (stroke_width / 2.0), n.y_coor).L(c2.x_coor + (stroke_width / 2.0), n.y_coor)
                genetree_svg.append(p)
                

        # We add the symbols
        for e, x_coor, y_coor in symbols:
            if e == "O":
                p = draw.Circle(cx=x_coor, cy=y_coor,r=stroke_width * symbol_scaling, stroke=stroke, stroke_width=0, fill='green')
                genetree_svg.append(p)
            if e == "D":
                p = draw.Circle(cx=x_coor, cy=y_coor,r=stroke_width * symbol_scaling, stroke=stroke, stroke_width=0, fill='yellow')
                genetree_svg.append(p)
            if e == "L":
                p = draw.Circle(cx=x_coor, cy=y_coor,r=stroke_width * symbol_scaling, stroke=stroke, stroke_width=0, fill='red')
                genetree_svg.append(p)
                
        for element in genetree_svg:            
            d.append(element)
        d.setRenderSize(h=self.height)
        d.saveSvg(url)
        return(d)
    
    def AnimateForward(self, url): 
        
        total_frames = 200
        
        d = self.d         
        p = draw.Rectangle(0,0, self.width, self.height, fill="black", id="Square")  
        d.append(p)        
                
        def draw_frame(i):            
            std.d.elements[-1].args["y"] += (self.height / total_frames)
            #d.setRenderSize(h=self.height)
            #d.saveSvg(url)
            return(d)    
        with draw.animate_jupyter(draw_frame, delay=0.05) as anim:
            for i in range(total_frames):
                anim.draw_frame(i)
    
           
    def ZoomOnPoint(self, point, total_frames, reset_view = False):

        d = self.d  

        x0,y0,x1,y1 = std.d.viewBox 
        
        x0f, y0f, x1f, y1f = point
        
        x0s = (x0f - x0) / float(total_frames)
        y0s = (y0f - y0) / float(total_frames)
        x1s = (x1f - x1) / float(total_frames)
        y1s = (y1f - y1) / float(total_frames)
        
                
        def draw_frame(i):
            std.d.viewBox = (x0 + i * x0s,
                             y0 + i * y0s,
                             x1 + i * x1s,
                             y1 + i * y1s,
                            )
                
            return(d)    
        with draw.animate_jupyter(draw_frame, delay=0.05) as anim:
            for i in range(total_frames):
                anim.draw_frame(i)
        
        if reset_view:
            std.d.viewBox = x0,y0,x1,y1
            
    def ZoomOnBranch(self, branch, total_frames, reset_view = False, margin=20):

        d = self.d  

        x0,y0,x1,y1 = std.d.viewBox 
        
        n = self.tree&branch
        
        x0f, y0f, x1f, y1f = n.x_coor - margin, n.y_coor - margin, n.x_coor + margin, n.y_coor + margin  
        
        x0s = (x0f - x0) / float(total_frames)
        y0s = (y0f - y0) / float(total_frames)
        x1s = (x1f - x1) / float(total_frames)
        y1s = (y1f - y1) / float(total_frames)
        
                
        def draw_frame(i):
            std.d.viewBox = (x0 + i * x0s,
                             y0 + i * y0s,
                             x1 + i * x1s,
                             y1 + i * y1s,
                            )
                
            return(d)    
        with draw.animate_jupyter(draw_frame, delay=0.05) as anim:
            for i in range(total_frames):
                anim.draw_frame(i)
        
        if reset_view:
            std.d.viewBox = x0,y0,x1,y1

    def draw_constraint(self, url, dn, rc, stroke="red", stroke_width=1):

        d = self.d          
        
        arrow = draw.Marker(-0.1, -0.5, 0.9, 0.5, scale=stroke_width*1.5, orient='auto')
        arrow.append(draw.Lines(-0.1, -0.5, -0.1, 0.5, 0.9, 0, fill=stroke, close=True))
        
        
        p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none', 
                                      stroke_dasharray="4", marker_end=arrow)  
            
        t = self.tree
        bdn = t&dn
        p_bdn = bdn.up
        
        brc = t&rc
        p_brc = brc.up
        
 
        if bdn.x_coor > brc.x_coor:    
 
            p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none', marker_end=arrow)  
            p.M(bdn.x_coor, bdn.y_coor).L(brc.x_coor, brc.y_coor)
            d.append(p)
        else:
 
            p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none', marker_end=arrow)  
            p.M(bdn.x_coor, bdn.y_coor).L(brc.x_coor, brc.y_coor)
            d.append(p)
        '''
        spacing = self.width / (self.n_leaves + 1)
        total_time = self.total_time
        width = self.width
        height = self.height
        d = self.d
        stroke_a = stroke_width / 2.0
       
        for n in self.tree.traverse():
            # Inner tree, vertical lines
            x_0, y_0, y_1 = n.x_coor, n.y_coor, self._y_scaling(n.dist)
            p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
            p.M(x_0, y_0).l(0, y_1)  
            d.append(p)                        
            if not n.is_leaf():        
                c1, c2 = n.get_children()                
                # Horizontal lines        
                # Inner Tree                
                p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none')  
                p.M(c1.x_coor - stroke_a, n.y_coor).L(c2.x_coor + stroke_a, n.y_coor)
                d.append(p)            
        d.setRenderSize(h=self.height)
        d.saveSvg(url)
        return(d)
    
    
        if node2events[n.name] == "T":                    
                        p = draw.Path(stroke=stroke, stroke_width=stroke_width, fill='none', 
                                      stroke_dasharray="4", marker_end=arrow)  

                        if c1.x_coor > c2.x_coor:                    
                            p.M(c1.x_coor, n.y_coor).L(c2.x_coor + radius, n.y_coor)
                        else:
                            p.M(c1.x_coor, n.y_coor).L(c2.x_coor - radius, n.y_coor)

    
        '''
        d.saveSvg(url)
        return(d)
            

    
    
            
