
warning('off', 'MATLAB:print:ContentTypeImageSuggested')

cd scripts

fig2
exportgraphics(fig, '../output/fig2.pdf', 'ContentType', 'vector')
    
fig3
exportgraphics(fig, '../output/fig3.pdf', 'ContentType', 'vector')

fig4
exportgraphics(fig, '../output/fig4.pdf', 'ContentType', 'vector')

fig5
exportgraphics(fig, '../output/fig5.pdf', 'ContentType', 'vector')

fig6
exportgraphics(fig, '../output/fig6.pdf', 'ContentType', 'vector')

cd ..

warning('on', 'MATLAB:print:ContentTypeImageSuggested')
