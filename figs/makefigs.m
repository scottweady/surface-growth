warning('off', 'MATLAB:print:ContentTypeImageSuggested')

cd scripts

for name = {'fig2', 'fig3', 'fig4', 'fig5', 'fig6'}
  run(name{1})
  exportgraphics(fig, fullfile('..', 'output', [name{1} '.pdf']), 'ContentType', 'vector')
  close
end

cd ..

warning('on', 'MATLAB:print:ContentTypeImageSuggested')
