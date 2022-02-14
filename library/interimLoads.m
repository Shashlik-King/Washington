function [loads] = interimLoads (plots,id,soil,loads,settings)
    if strcmp(id,'A04')
        if strcmp(settings.interface,'FAC') 
            loads.M=-654649;
            loads.H=10579;
            loads.Vc=29733;
            loads.Mz=21989;
        else
            loads.M=-489892;
            loads.H=7837;
            loads.Vc=21238;
            loads.Mz=15706;
        end			
   
    
   elseif strcmp(id,'B01')
        if strcmp(settings.interface,'FAC') 
            loads.M=-649435;
            loads.H=10533;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-486050;
            loads.H=7803;
            loads.Vc=1;
            loads.Mz=1;
        end
    elseif strcmp(id,'B02')
        if strcmp(settings.interface,'FAC') 
            loads.M=-650866;
            loads.H=10546;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-487104;
            loads.H=7812;
            loads.Vc=1;
            loads.Mz=1;
        end
   elseif strcmp(id,'B03')
        if strcmp(settings.interface,'FAC') 
            loads.M=-653627;
            loads.H=10570;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-489138;
            loads.H=7830;
            loads.Vc=1;
            loads.Mz=1;
        end          
 elseif strcmp(id,'B04')
        if strcmp(settings.interface,'FAC') 
            loads.M=-653831;
            loads.H=10572;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-489289;
            loads.H=7831;
            loads.Vc=1;
            loads.Mz=1;
        end    
 elseif strcmp(id,'C02')
        if strcmp(settings.interface,'FAC') 
            loads.M=-644261;
            loads.H=10485;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-482203;
            loads.H=7767;
            loads.Vc=1;
            loads.Mz=1;
        end 
 elseif strcmp(id,'C03')
        if strcmp(settings.interface,'FAC') 
            loads.M=-646776;
            loads.H=10509;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-484091;
            loads.H=7785;
            loads.Vc=1;
            loads.Mz=1;
        end 
  elseif strcmp(id,'C04')
        if strcmp(settings.interface,'FAC') 
            loads.M=-649435;
            loads.H=10533;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-486050;
            loads.H=7803;
            loads.Vc=1;
            loads.Mz=1;
        end 
  elseif strcmp(id,'D02')
        if strcmp(settings.interface,'FAC') 
            loads.M=-641342;
            loads.H=10457;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-480016;
            loads.H=7746;
            loads.Vc=1;
            loads.Mz=1;
        end 
 elseif strcmp(id,'D04')
        if strcmp(settings.interface,'FAC') 
            loads.M=-650457;
            loads.H=10542;
            loads.Vc=1;
            loads.Mz=1;
        else
            loads.M=-486803;
            loads.H=7809;
            loads.Vc=1;
            loads.Mz=1;
        end        
         
    
    
  elseif strcmp(id,'G01')
        if strcmp(settings.interface,'FAC')
            loads.M=-638715;
            loads.H=10431;
            loads.Vc=29214;
            loads.Mz=21989;
        else
            loads.M=-478048;
            loads.H=7727;
            loads.Vc=20867;
            loads.Mz=15706;
        end        
    elseif strcmp(id,'G04')
        if strcmp(settings.interface,'FAC')
            loads.M = -626358;
            loads.H = 10313;
            loads.Vc=28798;
            loads.Mz=21989;
        else
            loads.M = -468790;
            loads.H = 7639;
            loads.Vc=20570;
            loads.Mz=15706;
        end   
    elseif strcmp(id,'OSS')
        if strcmp(settings.interface,'FAC')
            loads.M = -1043474;
            loads.H = 22959;
            loads.Vc= 46824;
            loads.Mz= 20233;
        else
            loads.M = -829678;
            loads.H = 17658;
            loads.Vc= 38802;
            loads.Vt= 18802;
            loads.Mz= 14126;
        end   
    else
        if strcmp(settings.interface,'FAC') 
            loads.M=-545782;
            loads.H=7526;
            loads.Vc=29079;
            loads.Mz=21864;
        else
            loads.M=-409328;
            loads.H=5575;
            loads.Vc=20771;
            loads.Mz=15617;
        end	   
   end