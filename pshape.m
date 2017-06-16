figure
for i=1:N
plot(time, abs(p(i,:)))
hold on
end
head = '%d points';
head1 = sprintf(head, N);
title(head1)