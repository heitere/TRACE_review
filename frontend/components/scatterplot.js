'use client'

import { useCallback, useState, useEffect } from "react";
import { MyButton, ResetButton, AsyncButton } from "./buttons";
import Legend from "./legend";
import { scatterplot } from "./canvas";
import { EmbeddingScale, SettingsMenu } from "./PlotControls"
import SelectionInfo from "./selectionInfo";
import { Colorbar } from "./colorbar";
import CanvasWrapper from '@/components/canvas_wrapper'
import { ReactSelect } from "./utils"

// detault options
let pointSizeInitial = 5;
let maxNeighbors = 200;
const metrics = ["angular", "euclidean", "manhattan", "hamming", "dot"]
const metricOptions = metrics.map(v => { return { value: v, label: v } })
//const datasetOptions = ["Mouse CD45neg", "Mammoth", "GaussLine", "Human Immune", "Mouse Fibroblasts" ]
var datasetOptions = [];

let opacities;
let filteredPoints = [];
let numPoints;

const resetZoomHandler = () => {
    // scatterplot.reset();
    scatterplot.zoomToLocation([0.0, 0.0], 1.2,
        {
            transition: true,
            transitionDuration: 1000,
        });
}

const resetOpacityHandler = () => {
    scatterplot.set({
        opacityBy: 'density',
    })
}

const zoomToSelectionHandler = () => {
    const selected = scatterplot.get('selectedPoints')
    opacities.forEach((v, i) => { if (v == 1) { selected.push(i) } })
    console.log(`zooming to ${selected.length} points`)
    scatterplot.zoomToPoints(selected, {
        padding: 1.0,
        transition: true,
        transitionDuration: 1500,
    })
}

const printCurrentSelection = () => {
    let selectedPoints = scatterplot.get('selectedPoints')
    console.log(JSON.stringify(selectedPoints))
}


const showHDNeighbors = (embedding, pointColors, kNeighbors, metric) => {
    let selectedPoints = scatterplot.get('selectedPoints')

    return new Promise((resolve, reject) => {
        if (selectedPoints.length > 0) {
            fetch("/backend/computeHDNeighbors", {
                method: "POST",
                body: JSON.stringify({
                    k: kNeighbors,
                    points: selectedPoints,
                    hd_metric: metric,
                }),
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
                .then(res => res.json())
                .then(data => {
                    let zoomPoints = data["result"]
                    opacities = data["binary"]

                    scatterplot.draw(
                        {
                            x: embedding['x'],
                            y: embedding['y'],
                            z: pointColors["encoded_values"],
                            w: opacities,
                        },
                        {
                            // preventFilterReset: true // does not work!
                            filter: filteredPoints,
                            zDataType: pointColors["type"]
                        }
                    ).then(() => {
                        scatterplot.set({
                            opacityBy: 'w',
                            opacity: [0.03, 1],
                        })
                        resolve(true);
                        // if (zoomPoints.length > 0) {
                        //     scatterplot.zoomToPoints(zoomPoints, {
                        //         padding: 0.9,
                        //         transition: true,
                        //         transitionDuration: 1500,
                        //     })
                        // }
                    });
                    console.log(`showing HD-neighbors`)
                })
        } else {
            resolve(true);
        }
    })
}


const showHDReverseNeighbors = (embedding, pointColors, kNeighbors, metric) => {
    let selectedPoints = scatterplot.get('selectedPoints')

    return new Promise((resolve, reject) => {
        if (selectedPoints.length >= 1) {
            fetch("/backend/computeReverseHDNeighbors", {
                method: "POST",
                body: JSON.stringify({
                    k: kNeighbors,
                    points: selectedPoints,
                    hd_metric: metric
                }),
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
                .then(res => res.json())
                .then(data => {
                    opacities = data["binary"]

                    scatterplot.draw(
                        {
                            x: embedding['x'],
                            y: embedding['y'],
                            z: pointColors["encoded_values"],
                            w: opacities,
                        },
                        {
                            // preventFilterReset: true // does not work!
                            filter: filteredPoints,
                            zDataType: pointColors["type"]
                        }
                    ).then(() => {
                        scatterplot.set({
                            opacityBy: 'w',
                            opacity: [0.03, 1],
                        })
                        console.log(`showing reverse HD-neighbors`)
                        resolve(true);
                    })
                })
        } else {
            resolve(true);
        }
    })
}


const showIntrusions = (kNeighbors, metric) => {
    let selectedPoints = scatterplot.get('selectedPoints')

    return new Promise((resolve, reject) => {
        if (selectedPoints.length > 1) {
            fetch("/backend/intrusions", {
                method: "POST",
                body: JSON.stringify({
                    k: kNeighbors,
                    points: selectedPoints,
                    hd_metric: metric
                }),
                headers: {
                    "Content-type": "application/json; charset=UTF-8"
                }
            })
                .then(res => res.json())
                .then(data => {
                    let intrusions = data["result"]
                    scatterplot.deselect()
                    if (intrusions.length > 0) {
                        scatterplot.select(intrusions);
                        scatterplot.zoomToPoints(selectedPoints, {
                            padding: 0.9,
                            transition: true,
                            transitionDuration: 1500,
                        })
                    }
                    console.log(`showing ${intrusions.length} intrusions`)
                    resolve(true);
                })
        } else {
            resolve(true);
        }
    })
}

const filterPoints = (featureValues, filterValue) => {
    // only show points with value != 'filterValue'
    let show = [];
    filteredPoints.map((v) => {
        if (featureValues[v] != filterValue) {
            show.push(v)
        }
    });
    let hiddenLength = filteredPoints.length - show.length;
    filteredPoints = show;
    scatterplot.filter(show);
    console.log(`hiding ${filterValue} (${hiddenLength} points)`)
}

const unfilterPoints = (featureValues, filterValue) => {
    // show points with value 'filterValue'
    let show = [];
    featureValues.map((v, i) => {
        if (v == filterValue) show.push(i);
    })
    let addedLength = show.length;
    show = show.concat(filteredPoints)
    filteredPoints = show;
    scatterplot.filter(show)
    console.log(`showing ${filterValue} (${addedLength} points)`)
}

const resetPointFilter = () => {
    filteredPoints = [...Array(numPoints).keys()];
}

const getPointColors = (embName, featureName, embIndex) => {
    return new Promise((resolve, reject) => {
        fetch(`/backend/pointColor/${featureName}?embeddingName=${embName}&embeddingIndex=${embIndex}`)
            .then(res => res.json())
            .then(res => {
                resolve(res)
            })
    })
}

export const isSameElements = (a, b) => {
    if (a.length !== b.length) return false;
    const aSet = new Set(a);
    return b.every((value) => aSet.has(value));
};

const showEmbedding = (embedding, pointColor, useTransition = true, preventFilterReset = true) => {
    let cview = scatterplot.get('cameraView');
    if (!preventFilterReset) resetPointFilter();

    scatterplot
        .draw(
            {
                x: embedding['x'],
                y: embedding['y'],
                z: pointColor["encoded_values"],
                w: opacities
            },
            {
                transition: useTransition,
                transitionDuration: 1500,
                transitionEasing: 'quadInOut',
                // preventFilterReset: preventFilterReset,
                filter: preventFilterReset ? filteredPoints : undefined,
                zDataType: pointColor["type"]
            }
        ).then(() => {
            scatterplot.set({
                cameraView: cview,
                pointColor: Object.values(pointColor["colorMap"]),
                pointColorActive: "#ffba08", //"#55308d"
            })
        })
    scatterplot.set({
        opacityBy: scatterplot.get('opacityBy'),
    })
}


const fetchEmbedding = (embName, embScale) => {
    return new Promise((resolve, reject) => {
        fetch(`/backend/embedding?embName=${embName}&embScale=${embScale}`)
            .then(res => res.json())
            .then(res => {
                resolve(res)
            })
    })
}


export default function Scatterplot() {
    const [isLoading, setIsLoading] = useState(true);
    const [isLoadingData, setIsLoadingData] = useState(false);

    const [embeddingName, setEmbeddingName] = useState("UMAP");
    const [activeEmbedding, setActiveEmbedding] = useState(null)
    const [embeddingOptions, setEmbeddingOptions] = useState(null);
    const [selectedScale, setScale] = useState({ "UMAP": 0 });

    const [legendVisibility, setLegendVisibility] = useState("visible");
    const [kNeighbors, setkNeighbors] = useState(50);
    const [pointColorOptions, setPointColorOptions] = useState(null);
    const [selectedPointColor, setSelectedPointColor] = useState("initialPointColors");
    const [pointColors, setPointColors] = useState(
        {
            "values": 0,
            "encoded_values": 0,
            "colorMap": { "none": "#444444" },
            "type": "categorical"
        }
    );
    const [scatterLoaded, setScatterLoaded] = useState(false);
    const [scatterplotState, setScatterplot] = useState(null);
    const [datasetInfo, setDatasetInfo] = useState("");
    const [selectedMetric, setMetric] = useState("angular");
    const [datasetName, setDatasetName] = useState(datasetOptions[0]);

    const handleDatasetSelect = (newDatasetName) => {
        console.log(`handleDatasetSelect ${newDatasetName}`)
        setIsLoadingData(true);
        setDatasetName(newDatasetName);
    }

    const loadDatasetConfiguration = (datasetName, unsetLoadingFn) => {
        return new Promise((resolve, reject) => {
            let currScale = 0;
            var newEmbeddingName;

            fetch(`/backend/loadDataset?datasetName=${datasetName}`)
                .then(res => res.json())
                .then(res => {
                    setMetric(res["hd_metric"]);
                    setDatasetInfo(res["dataset_info"]);
                    setDatasetName(res["dataset_name"]);

                    // fetch options and scales
                    fetch("/backend/embeddingOptions")
                        .then(res => res.json())
                        .then(res => {
                            setEmbeddingOptions(res)
                            newEmbeddingName = Object.keys(res)[0];
                            setEmbeddingName(newEmbeddingName);

                            const scaleObj = {};
                            for (const key of Object.keys(res)) {
                                scaleObj[key] = 0;
                            }
                            setScale(scaleObj);

                            Promise.all([fetchEmbedding(newEmbeddingName, currScale), fetch("/backend/pointColorOptions")])
                                .then(([newEmbedding, res]) => {
                                    setActiveEmbedding(newEmbedding)
                                    numPoints = newEmbedding["x"].length;
                                    filteredPoints = [...Array(numPoints).keys()];
                                    opacities = new Array(numPoints).fill(0);
                                    console.log(`number of points ${numPoints}`)

                                    res.json().then(res => {
                                        let newSelectedPointColor = res["result"][0]
                                        let newPointColorOptions = ["none"].concat(res["result"])
                                        setPointColorOptions(newPointColorOptions.map(v => { return { value: v, label: v } }));
                                        setSelectedPointColor(newSelectedPointColor);
                                        console.log(`set pointColors to ${newSelectedPointColor}`)

                                        getPointColors(newEmbeddingName, newSelectedPointColor, currScale)
                                            .then((pc) => {
                                                setPointColors(pc);
                                                console.log("Finished loading data")
                                                unsetLoadingFn(false);
                                                resolve(true);
                                            })
                                    })
                                })
                        })
                })
                .catch(err => {
                    console.log(err)
                    reject(err)
                }
                );
        })
    }


    const handlePointColorSelect = (newPointColor) => {
        console.log(`handlePointColorSelect ${newPointColor}`)
        setSelectedPointColor(newPointColor);
        getPointColors(embeddingName, newPointColor, selectedScale[embeddingName])
            .then((res) => {
                setPointColors(res);
                showEmbedding(activeEmbedding, res, false, false);
            })
    }

    // fetch the list of unstable points from the backend
    const showUnstablePoints = () => {
        return new Promise((resolve, reject) => {
            scatterplot.deselect();
            fetch(`/backend/getUnstablePoints?embeddingName=${embeddingName}&embeddingScale=${selectedScale[embeddingName]}`)
                .then(res => res.json())
                .then(res => {
                    console.log(`showing ${res['result'].length} unstable points`)
                    scatterplot.select(res['result']);
                    resolve(true);
                })
        })
    }



    const handleEmbeddingSelect = (newEmbeddingName) => {
        setEmbeddingName(newEmbeddingName);
        fetchEmbedding(newEmbeddingName, selectedScale[newEmbeddingName])
            .then(newEmbedding => {
                setActiveEmbedding(newEmbedding);

                if (selectedPointColor.includes('quality')) {
                    // pointColor was quality of old embedding ... recompute
                    getPointColors(newEmbeddingName, selectedPointColor, selectedScale[newEmbeddingName])
                        .then((res) => {
                            if ("none" in res["colorMap"]) {
                                setSelectedPointColor("none");
                            }
                            setPointColors(res);
                            showEmbedding(newEmbedding, res);
                        })
                } else {
                    showEmbedding(newEmbedding, pointColors);
                }
            })
    }

    const handleScaleSelect = (scale) => {
        setScale(prevScale => {
            let newScale = { ...prevScale }
            newScale[embeddingName] = scale
            return newScale
        });
        fetchEmbedding(embeddingName, scale)
            .then(newEmbedding => {
                setActiveEmbedding(newEmbedding);
                if (selectedPointColor.includes('quality')) {
                    // pointColor was quality of old embedding ... recompute
                    console.log("recomputing point colors")
                    getPointColors(embeddingName, selectedPointColor, scale)
                        .then((res) => {
                            if ("none" in res["colorMap"]) {
                                setSelectedPointColor("none");
                            }
                            setPointColors(res);
                            showEmbedding(newEmbedding, res);
                        })
                } else {
                    showEmbedding(newEmbedding, pointColors);
                }
            })
    }

    useEffect(() => {
        async function fetchData() {
            try {
                loadDatasetConfiguration(datasetOptions[0], setIsLoading);
            } catch (error) {
                console.error('Error fetching data:', error);
            }
        }

        // fetch options and scales
        fetch("/backend/datasetOptions")
            .then(res => res.json())
            .then(res => {
                datasetOptions = res["result"]
                setDatasetName(datasetOptions[0])
                fetchData();
            })
    }, []);


    useEffect(() => {
        if (isLoadingData) {
            loadDatasetConfiguration(datasetName, setIsLoadingData)
        }
    }, [isLoadingData]);

    useEffect(() => {
        if (scatterLoaded && !isLoading && !isLoadingData) {
            console.log(`drawing initial embedding ${embeddingName}`)
            showEmbedding(activeEmbedding, pointColors, false, false);
            scatterplot.deselect();
            resetOpacityHandler();
        }
    }, [scatterLoaded, isLoading, isLoadingData]);

    // need to handle case when data changes, i.e. point color options, .... everything must be updated
    if (isLoading) {
        return (
            <div role="status" className="absolute w-screen h-screen bg-white bg-opacity-75 flex justify-center items-center">
                <svg aria-hidden="true" className="w-8 h-8 text-gray-200 animate-spin dark:text-gray-600 fill-blue-600" viewBox="0 0 100 101" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M100 50.5908C100 78.2051 77.6142 100.591 50 100.591C22.3858 100.591 0 78.2051 0 50.5908C0 22.9766 22.3858 0.59082 50 0.59082C77.6142 0.59082 100 22.9766 100 50.5908ZM9.08144 50.5908C9.08144 73.1895 27.4013 91.5094 50 91.5094C72.5987 91.5094 90.9186 73.1895 90.9186 50.5908C90.9186 27.9921 72.5987 9.67226 50 9.67226C27.4013 9.67226 9.08144 27.9921 9.08144 50.5908Z" fill="currentColor" /><path d="M93.9676 39.0409C96.393 38.4038 97.8624 35.9116 97.0079 33.5539C95.2932 28.8227 92.871 24.3692 89.8167 20.348C85.8452 15.1192 80.8826 10.7238 75.2124 7.41289C69.5422 4.10194 63.2754 1.94025 56.7698 1.05124C51.7666 0.367541 46.6976 0.446843 41.7345 1.27873C39.2613 1.69328 37.813 4.19778 38.4501 6.62326C39.0873 9.04874 41.5694 10.4717 44.0505 10.1071C47.8511 9.54855 51.7191 9.52689 55.5402 10.0491C60.8642 10.7766 65.9928 12.5457 70.6331 15.2552C75.2735 17.9648 79.3347 21.5619 82.5849 25.841C84.9175 28.9121 86.7997 32.2913 88.1811 35.8758C89.083 38.2158 91.5421 39.6781 93.9676 39.0409Z" fill="currentFill" /></svg>
                <span className="sr-only">Loading...</span>
            </div>
        );
    }

    return (
        <>
            <div className="flex-grow h-screen max-h-screen bg-white relative min-w-[400px]">
                <CanvasWrapper setScatterLoaded={setScatterLoaded} setScatterplot={setScatterplot} />
                <div className="absolute flex flex-wrap top-0 my-1 mr-6 ml- items-center display-block">
                    <MyButton onClick={() => resetOpacityHandler()}>reset opacity</MyButton>
                    <AsyncButton onClick={
                        () => showHDNeighbors(activeEmbedding, pointColors, kNeighbors, selectedMetric)
                    }>HD neighbors</AsyncButton>
                    {/*<AsyncButton onClick={
                        () => showHDReverseNeighbors(activeEmbedding, pointColors, kNeighbors, selectedMetric)
                    }>reverse neighbors</AsyncButton>*/}
                    <AsyncButton onClick={() => showIntrusions(kNeighbors, selectedMetric)}>intrusions</AsyncButton>
                    <ResetButton onClick={resetZoomHandler} />
                </div>
            </div>
            <SettingsMenu
                scatterplot={scatterplotState}
                selectedPointColor={selectedPointColor}
                pointColors={pointColors}
                pointSizeInitial={pointSizeInitial}
                pointColorOptions={pointColorOptions}
                pointColorOnChange={handlePointColorSelect}
                setLegendVisibility={setLegendVisibility}
                legendVisibility={legendVisibility}
                kNeighbors={kNeighbors}
                maxNeighbors={maxNeighbors}
                handlekNeighborSelect={setkNeighbors}
                metricOptions={metricOptions}
                selectedMetric={selectedMetric}
                metricOnChange={setMetric}
                showUnstablePoints={showUnstablePoints}
            >
                {/* Dataset */}
                <div className="flex flex-wrap items-center p-3 justify-between">
                    <label className="text-sm text-gray-500 w-fit min-w-fit mr-2" htmlFor='pointColorSelect'>Data</label>
                    <ReactSelect
                        options={datasetOptions.map(v => { return { value: v, label: v } })}
                        selected={datasetName}
                        onChange={handleDatasetSelect} />
                </div>


                {/* Embedding method */}
                <div className="flex flex-wrap items-center p-3 justify-between">
                    <label className="text-sm text-gray-500 w-fit min-w-fit mr-2" htmlFor='pointColorSelect'>DR method</label>
                    <ReactSelect
                        options={Object.keys(embeddingOptions).map(v => { return { value: v, label: v } })}
                        selected={embeddingName}
                        onChange={handleEmbeddingSelect} />
                </div>

                {/* Embedding scale (neighborhood size) */}
                <div className="flex flex-wrap items-center p-3 justify-between">
                    <label className="text-sm text-gray-500 w-fit min-w-fit mr-2">emb scale</label>
                    <EmbeddingScale
                        // how many embeddings are there for this method?
                        maxValue={embeddingOptions[embeddingName] - 1}
                        currValue={selectedScale[embeddingName]}
                        onChange={handleScaleSelect}
                    />
                </div>

            </SettingsMenu>
            <div className="fixed bottom-0 left-0 flex-wrap flex-row m-2 max-h-[90%] overflow-auto">
                {
                    pointColors["type"] === "categorical" ?
                        <Legend
                            colormap={pointColors["colorMap"]}
                            title={selectedPointColor}
                            filter={(filterValue) => filterPoints(pointColors["values"], filterValue)}
                            unfilter={(filterValue) => unfilterPoints(pointColors["values"], filterValue)}
                            visibility={legendVisibility}
                        /> :
                        <Colorbar
                            colormap={pointColors["colorMap"]}
                            title={selectedPointColor}
                            visibility={legendVisibility}
                        />
                }
            </div>

            <SelectionInfo scatterplot={scatterplotState} numPoints={numPoints} datasetInfo={datasetInfo} />

            <div
                //visibility={isLoadingData ? "visible" : "hidden"}
                role="status" className="absolute w-screen h-screen bg-white bg-opacity-75 flex justify-center items-center"
                style={{ visibility: isLoadingData ? "visible" : "hidden" }}
            >
                <svg aria-hidden="true" className="w-8 h-8 text-gray-200 animate-spin dark:text-gray-600 fill-blue-600" viewBox="0 0 100 101" fill="none" xmlns="http://www.w3.org/2000/svg"><path d="M100 50.5908C100 78.2051 77.6142 100.591 50 100.591C22.3858 100.591 0 78.2051 0 50.5908C0 22.9766 22.3858 0.59082 50 0.59082C77.6142 0.59082 100 22.9766 100 50.5908ZM9.08144 50.5908C9.08144 73.1895 27.4013 91.5094 50 91.5094C72.5987 91.5094 90.9186 73.1895 90.9186 50.5908C90.9186 27.9921 72.5987 9.67226 50 9.67226C27.4013 9.67226 9.08144 27.9921 9.08144 50.5908Z" fill="currentColor" /><path d="M93.9676 39.0409C96.393 38.4038 97.8624 35.9116 97.0079 33.5539C95.2932 28.8227 92.871 24.3692 89.8167 20.348C85.8452 15.1192 80.8826 10.7238 75.2124 7.41289C69.5422 4.10194 63.2754 1.94025 56.7698 1.05124C51.7666 0.367541 46.6976 0.446843 41.7345 1.27873C39.2613 1.69328 37.813 4.19778 38.4501 6.62326C39.0873 9.04874 41.5694 10.4717 44.0505 10.1071C47.8511 9.54855 51.7191 9.52689 55.5402 10.0491C60.8642 10.7766 65.9928 12.5457 70.6331 15.2552C75.2735 17.9648 79.3347 21.5619 82.5849 25.841C84.9175 28.9121 86.7997 32.2913 88.1811 35.8758C89.083 38.2158 91.5421 39.6781 93.9676 39.0409Z" fill="currentFill" /></svg>
                <span className="sr-only">Loading...</span>
            </div>
        </>
    );
}


