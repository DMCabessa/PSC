cd data
load horse-colic(processed).test ;

for i = 1:size(horse_colic_processed_,2)
   col = horse_colic_processed_((horse_colic_processed_(:,i) ~= -1),i) ;
   avg = mean(col) ;
   for j = 1:size(horse_colic_processed_,1)
       if horse_colic_processed_(j,i) == -1
           horse_colic_processed_(j,i) = avg ;
       end % if horse_colic_processed_
   end % for j
end % for i

fid = fopen('horse-colic(post-processed).test','w');

for i = 1:size(horse_colic_processed_,1)
	for j = 1:size(horse_colic_processed_,2)
		fprintf(fid,'%.2f, ',horse_colic_processed_(i,j)) ;
	end % for j
	fprintf(fid,'\n') ;
end % for i

cd ..