classdef cls_machine
    %UNTITLED10 此处提供此类的摘要
    %   此处提供详细说明

    properties
        l_stk
        r_outgap
        num_slot
        num_poles
        ang_rad_toothsurf
        num_fea_step
        slot_erea
        num_ele_layer
        num_nd_layer
    end

    methods
        function obj = cls_machine()
            obj.l_stk = 95;
            obj.r_outgap = 49;
            obj.num_slot = 12;
            obj.num_poles = 10;  
            obj.num_fea_step = 120;  % number of steps in EM FEA
            obj.num_ele_layer = 40;
            obj.num_nd_layer = obj.num_ele_layer+1;
            obj.ang_rad_toothsurf = [16.3448 43.6552]/180*pi;  % the angular position of two sides of the first tooth
            obj.slot_erea=zeros(obj.num_slot,2);
            tooth_ang = 2*pi/obj.num_slot;
            for i=1:obj.num_slot
                if i==1
                    obj.slot_erea(i,:)= [obj.ang_rad_toothsurf(2), obj.ang_rad_toothsurf(1)+tooth_ang];
                else
                    obj.slot_erea(i,:)= obj.slot_erea(1,:)+tooth_ang*(i-1);
                end
            end
            obj.slot_erea = mod(obj.slot_erea,2*pi);
        end
    end
end