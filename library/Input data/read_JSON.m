function [output] = read_JSON(fname)
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    output = jsondecode(str);
end