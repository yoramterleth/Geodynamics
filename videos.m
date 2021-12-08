%% makin g vieos 

% give v1 an array with the frames ( make using, in a for loop:
% F = getframe ;
% Frame_torque(i)= F;
% ) 
% vname should be the name of the outfile (and optional file path) 

function [] = videos(v1, vname)


%% make video 1 

% create the video writer with 1 fps
  writerObj = VideoWriter([vname,'.avi']);
  writerObj.FrameRate = .5;
  % set the seconds per image
%   if isempty(writerObj.cdata)
%       writerObj.cdata =  
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(v1)
    % convert the image to a frame
    frame = v1(i) ;    
    writeVideo(writerObj, frame);
end
%close the writer object
close(writerObj);

end 