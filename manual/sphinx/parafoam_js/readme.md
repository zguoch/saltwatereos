**Using [vtk.js](https://kitware.github.io/vtk-js/) present CFD results on web pages.**

# Generic js library development: view vtkjs scene

The basic develop procedure can be found on [vtk.js Doc](https://kitware.github.io/vtk-js/docs/intro_vtk_as_external_script.html).

Here I write a js file `src/index.js` as the main interface to read **scene** file and render it on web page.
Code is simple ans shown below.

```js
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
    // --------------------------------------------------------------------------------------------
    // load data and display. Actually .vtkjs is a zip file just with extentation name of vtkjs
    // --------------------------------------------------------------------------------------------
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
const userParams = vtkURLExtract.extractURLParameters();
if (userParams.url || userParams.fileURL) {
    showScene(userParams);
}
```

## Generate target library js file and test

- `npm run build`
- `npm start`

The online **scene explorer** can be found [here](https://kitware.github.io/vtk-js/examples/SceneExplorer.html)

Note that the file name must use https port, if you double click the index.html file, the .vtkjs file can not be loaded.

# Usage in sphinx manual

## File path

- _static/vtk_js/parafoam.js
- _static/vtk_js/index.html
- _static/vtk_js/xxx.vtkjs

## Use in rst file

The file name is specified through url, using `?` follow `_static/vtk_js/index.html`, and filename `fileURL=xxx.vtkjs`  follow `?`. see example: 

```rst
`3D model <../_static/vtk_js/index.html?fileURL=3D_TwoLayers.vtkjs>`_
```




