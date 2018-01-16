function [  ] = show_surf( grid_size,func,name,video_name,max_val )
%SHOW_surf get the demo video

x = 0:1/grid_size:1;
y = 0:1/grid_size:1;

[X,Y] = meshgrid(x,y);



time_step=0.001;

aviobj=VideoWriter(video_name);
aviobj.open();

for i = 1:100
    answer = func(i*time_step,grid_size);
    surf(X,Y,answer);
    title(name);
    axis([0 1,0 1,-max_val max_val]);
    caxis([0,max_val]);
    shading interp
    drawnow;
    frame=getframe(gcf);
    writeVideo(aviobj,frame);
    pause(time_step);
end

close(aviobj);

end

