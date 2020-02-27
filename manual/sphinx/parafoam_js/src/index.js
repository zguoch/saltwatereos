import 'vtk.js/Sources/favicon';

import vtkFullScreenRenderWindow from 'vtk.js/Sources/Rendering/Misc/FullScreenRenderWindow';
import vtkHttpSceneLoader from 'vtk.js/Sources/IO/Core/HttpSceneLoader';
import DataAccessHelper from 'vtk.js/Sources/IO/Core/DataAccessHelper';
import HttpDataAccessHelper from 'vtk.js/Sources/IO/Core/DataAccessHelper/HttpDataAccessHelper';
import vtkURLExtract from 'vtk.js/Sources/Common/Core/URLExtract';

function showScene(option) {
    // ----------------------------------------------------------------------------
    // Standard rendering code setup
    // ----------------------------------------------------------------------------

    const fullScreenRenderer = vtkFullScreenRenderWindow.newInstance();
    const renderer = fullScreenRenderer.getRenderer();
    const renderWindow = fullScreenRenderer.getRenderWindow();

    // ----------------------------------------------------------------------------
    // load data and display. Actually .vtkjs is a zip file 
    // just with extentation name of vtkjs
    // ----------------------------------------------------------------------------
    HttpDataAccessHelper.fetchBinary(option.fileURL).then((zipContent) => {
        // container.removeChild(progressContainer);
        const dataAccessHelper = DataAccessHelper.get('zip', {
            zipContent,
            callback: (zip) => {
                const sceneImporter = vtkHttpSceneLoader.newInstance({
                    renderer,
                    dataAccessHelper,
                });
                sceneImporter.setUrl('index.json');
                sceneImporter.onReady(() => {
                    renderWindow.render();
                });
            },
        });
    });
}

// showScene();

const userParams = vtkURLExtract.extractURLParameters();
if (userParams.url || userParams.fileURL) {
    showScene(userParams);
}

