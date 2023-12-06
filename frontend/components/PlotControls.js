import { useState} from "react";
import NamedSlider, { RadioSelect, ReactSelect } from "./utils";
import { MyButton, SettingsButton, ChevronRightButton} from "./buttons";
import { Switch } from '@headlessui/react'
import { ChevronDownIcon, ChevronUpIcon } from "@heroicons/react/24/solid"
import { Infobox } from "./infobox";


export function EmbeddingSelection(props) {
  const { onChange, defaultValue, options, ...other } = props;
  const [visibility, setVisibility] = useState("hidden")

  const toggleVisibility = () => {
    if (visibility == "visible") {
      setVisibility("hidden")
    }
    else {
      setVisibility("visible")
    }
  }

  return (
    <div className="w-fit h-fit mr-2">
      <MyButton onClick={toggleVisibility}>
        {<>
          <p className="select-none">DR method</p>
          {visibility == "visible" ? (<ChevronUpIcon className="w-4 h-4 ml-2" />) : (<ChevronDownIcon className="w-4 h-4 ml-2" />)}
        </>
        }
      </MyButton>
      <div className={visibility + " z-50 text-sm text-gray-500 absolute flex flex-col select-none w-fit bg-gray-100/90 rounded-lg bg-clip-border mt-1 ml-1  px-5 pb-2 pt-2.5"}>
        <RadioSelect onChange={onChange} options={options} defaultValue={defaultValue} />
      </div>
    </div>
  )
}

export function EmbeddingScale(props) {
  const {
    onChange,
    maxValue,
    currValue
  } = props

  return (
    <div className="w-1/2 min-w-[100px]">
      <div className="flex justify-between text-xs text-gray-500">
        <span className="inline select-none">global</span>
        <span className="inline select-none">local</span>
      </div>
      {
        maxValue < 1 ?
          <input
            className="transparent h-[2px] w-full cursor-pointer appearance-none border-transparent bg-neutral-300 my-2"
            type="range" min={0} max={maxValue} step={1} value={currValue}
            onChange={(event) => onChange(+event.target.value)} disabled
          />
          :
          <input
            className="transparent h-[2px] w-full cursor-pointer appearance-none border-transparent bg-neutral-300 my-2"
            type="range" min={0} max={maxValue} step={1} value={currValue}
            onChange={(event) => onChange(+event.target.value)}
          />
      }
    </div>
  )
}

const computeNeighbors = (kNeighbors, metric) => {
  return new Promise((resolve, reject) => {
    fetch(`/backend/precomputeAllNeighbors?maxK=${kNeighbors}&hd_metric=${metric}`, {
      method: "POST",
      body: {},
      headers: {
        "Content-type": "application/json; charset=UTF-8"
      }
    })
      .then(res => res.json())
      .then(data => {
        resolve(true);
      })
  })
}


export function SettingsMenu(props) {
  const {
    scatterplot,
    selectedPointColor,
    pointColors,
    pointSizeInitial,
    pointColorOptions,
    pointColorOnChange,
    setLegendVisibility,
    legendVisibility,
    kNeighbors,
    maxNeighbors,
    handlekNeighborSelect,
    metricOptions,
    selectedMetric,
    metricOnChange,
    showUnstablePoints,
    children } = props

  const [visibility, setVisibility] = useState('visible')
  const [selectedPointSize, setPointSize] = useState(pointSizeInitial);

  const toggleVisibility = () => {
    if (visibility == "visible") setVisibility("hidden"); else setVisibility("visible");
  }

  const handlePointSizeSelect = (pointSize) => {
    setPointSize(pointSize);
    scatterplot.set({ pointSize });
  };

  const toggleLegendVisibility = () => {
    if (legendVisibility == "visible") {
      setLegendVisibility("hidden")
    } else {
      setLegendVisibility("visible")
    }
  }

  const handleOpacitySelect = (opacity) => {
    console.log(`opacity by is ${scatterplot.get('opacityBy')}`)
    if (scatterplot.get('opacityBy') == 'valueW') {
      scatterplot.set({
        opacityBy: 'w',
        opacity: [opacity, 1],
      })
    } else {
      scatterplot.set({
        opacityBy: null,
        opacity: opacity,
      })
    }
  }


  if (visibility == "visible") {
    return (
      <>
        <div className="select-none right-0 w-1/4 min-w-[200px] h-screen bg-gray-100">
          <span className="relative">
            <div className="absolute -left-[20px] top-2">
              <ChevronRightButton onClick={toggleVisibility} />
            </div>
          </span>

          <div className="flex flex-col min-h-fit items-top justify-center text-center p-5 overflow-auto">

            <h3 className="text-lg font-medium leading-6 text-gray-900 w-fit mt-8" >
              Settings
            </h3>

            {children}

            {/* Point Colors */}
            {console.log(`selected point color is ${selectedPointColor}`)}
            <div className="flex flex-wrap items-center p-3 justify-between">
              <label className="text-sm text-gray-500 w-fit min-w-fit mr-2" htmlFor='pointColorSelect'>point colors</label>
              <ReactSelect
                options={pointColorOptions}
                selected={selectedPointColor}
                onChange={pointColorOnChange} />
            </div>

            {/* Point Size */}
            <div className='flex flex-wrap items-center p-3 justify-between'>
              <label className="text-sm text-gray-500 w-fit min-w-fit mr-2">point size</label>
              <input
                className="transparent h-[3px] w-1/2 min-w-[100px] cursor-pointer appearance-none border-transparent bg-neutral-300 my-2"
                type="range" min={1} max={10} step={1} defaultValue={selectedPointSize}
                onChange={(event) => handlePointSizeSelect(+event.target.value)} />
            </div>

            {/* Point Opacity */}
            <div className='flex flex-wrap items-center p-3 justify-between'>
              <label className="text-sm text-gray-500 w-fit min-w-fit mr-2">point opacity</label>
              <input
                className="transparent h-[3px] w-1/2 min-w-[100px] cursor-pointer appearance-none border-transparent bg-neutral-300 my-2"
                type="range" min={0} max={0.6} step={0.005} defaultValue={0.4}
                onChange={(event) => handleOpacitySelect(+event.target.value)} />
            </div>

            {/* Infobox */}
            <Infobox scatterplot={scatterplot} selectedPointColor={selectedPointColor} pointColors={pointColors} />

            {/* Legend */}
            <Switch.Group>
              <div className="flex flex-wrap items-center p-3 justify-between">
                <Switch.Label className="text-sm text-gray-500 w-fit min-w-fit mr-2" htmlFor='hoverSwitch'>show legend</Switch.Label>
                <div className="w-1/2 flex items-start justify-start">
                  <Switch
                    id='hoverSwitch'
                    checked={legendVisibility == "visible" ? true : false}
                    onChange={toggleLegendVisibility}
                    className={`${legendVisibility == "visible" ? 'bg-blue-600' : 'bg-gray-200'
                      } relative inline-flex h-6 w-11 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2`}
                  >
                    <span
                      className={`${legendVisibility == "visible" ? 'translate-x-6' : 'translate-x-1'
                        } inline-block h-4 w-4 transform rounded-full bg-white transition-transform`}
                    />
                  </Switch>
                </div>
              </div>
            </Switch.Group>


            {/* Neighbors */}
            <h3 className="text-lg font-medium leading-6 text-gray-900 w-fit mt-8" >
              Embedding Quality
            </h3>

            {/* Distance measures */}
            <div className="flex flex-wrap items-center p-3 justify-between">
              <label className="text-sm text-gray-500 w-fit min-w-fit mr-2" htmlFor='pointColorSelect'>HD metric</label>
              <ReactSelect options={metricOptions} selected={selectedMetric} onChange={metricOnChange}
                menuPlacement={'top'} />
            </div>

            {/* kNN */}
            <div className='flex flex-wrap items-center p-3 justify-between'>
              <NamedSlider onChange={handlekNeighborSelect} defaultValue={kNeighbors} min="0" max={maxNeighbors} step="10" label="neighbors" />
            </div>

            {/* Compute Neighbors */}
            {/* <AsyncButton onClick={() => computeNeighbors(maxNeighbors, selectedMetric)}>precompute neighbors</AsyncButton> */}
            {/* <AsyncButton onClick={() => showUnstablePoints()}>show unstable points</AsyncButton> */}
          </div>
        </div>
      </>
    )
  } else {
    return (
      <>
        <div className="fixed right-0 top-0 my-2">
          <SettingsButton onClick={toggleVisibility} />
        </div>
      </>
    )
  }
}