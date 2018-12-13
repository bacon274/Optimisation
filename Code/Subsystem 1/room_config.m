classdef room_config
    properties
        window_1 = window(1.44,1); % note: all u values set to 1 right now, could be much higher
        window_2 = window(1.44,1);
        window_3 = window(0.15,1);
        window_4 = window(0.15,1);
        window_5 = window(1.05,1);
        Panel_data_power = [104,88,62,42,0];
        Panel_data_g = [0,0.24,0.31,0.38,0.55];
        Panel_data_c = [600,600,600,600,0];
        Panel_transparency_list = [0,10,30,50,100] 
        windows
    end
    methods
        function obj = room_config(C)
            
            obj.window_1.g = obj.Panel_data_g(C(1)+1);
            obj.window_1.power = obj.Panel_data_power(C(1)+1);
            obj.window_1.cost = obj.Panel_data_c(C(1)+1);
            obj.window_1.transparency = obj.Panel_transparency_list(C(1)+1);
            
            obj.window_2.g = obj.Panel_data_g(C(2)+1);
            obj.window_2.power = obj.Panel_data_power(C(2)+1);
            obj.window_2.cost = obj.Panel_data_c(C(2)+1);
            obj.window_2.transparency = obj.Panel_transparency_list(C(2)+1);
            
            obj.window_3.g = obj.Panel_data_g(C(3)+1);
            obj.window_3.power = obj.Panel_data_power(C(3)+1);
            obj.window_3.cost = obj.Panel_data_c(C(3)+1);
            obj.window_3.transparency = obj.Panel_transparency_list(C(3)+1);
            
            obj.window_4.g = obj.Panel_data_g(C(4)+1);
            obj.window_4.power = obj.Panel_data_power(C(4)+1);
            obj.window_4.cost = obj.Panel_data_c(C(4)+1);
            obj.window_4.transparency = obj.Panel_transparency_list(C(4)+1);
            
            obj.window_5.g = obj.Panel_data_g(C(5)+1);
            obj.window_5.power = obj.Panel_data_power(C(5)+1);
            obj.window_5.cost = obj.Panel_data_c(C(5)+1);
            obj.window_5.transparency = obj.Panel_transparency_list(C(5)+1);
            
            obj.windows = C(:);
           
           
                                       %C is a combination of window values 0 to 4, corresponding to g values
     
            
        end
    end
end


%{

            obj.window_1.g = C(1)/10;
            obj.window_2.g = C(2)/10;
            obj.window_3.g = C(3)/10;
            obj.window_4.g = C(4)/10;
            obj.window_5.g = C(5)/10;
%}