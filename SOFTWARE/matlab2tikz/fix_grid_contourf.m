function fix_grid_contourf(filepath)

    if(exist(filepath,'file') == false)
        error('File does not exist')
    end
    
    fid = fopen(filepath,'r');
    C = textscan(fid,'%s','delimiter', '\n');
    fclose(fid);
    
    
    begin_axis = find(strcmp(C{1,:},'\begin{axis}[%')) ;
    end_axis = min(find(strcmp(C{1,:},']')));
    
    loc = find(strcmp(C{1,:},'\end{tikzpicture}%')) - 1;
    
    C = C{1};
    % write data back to file:
    fid = fopen(filepath,'w');
    fprintf(fid,'%s \n',C{1:loc});
    fprintf(fid,'%s \n',C{begin_axis:end_axis});
    fprintf(fid,'%s \n','\end{axis}');
    fprintf(fid,'%s',C{loc+1});
    fclose(fid);
    
end
