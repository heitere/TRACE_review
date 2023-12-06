import { useState, useEffect } from "react";
import { InformationCircleIcon } from "@heroicons/react/20/solid"


// https://tailwindcomponents.com/component/profile-information-card-horizon-ui-tailwind
export default function SelectionInfo(props) {
    const { scatterplot, numPoints, datasetInfo } = props;

    const [numSelected, setNumSelected] = useState(0);
    useEffect(() => {
        if (scatterplot != null) {
            console.log("subscribing to select action");
            scatterplot.subscribe('select', onSelect);
            scatterplot.subscribe('deselect', onDeselect);
        } else {
            console.log("Info: scatterplot is null");
        }
    }, [scatterplot]);

    const onSelect = (points) => {
        if (points["points"].length > 0) {
            setNumSelected(points["points"].length)
        }
    }

    const onDeselect = () => {
        let cview = scatterplot.get('cameraView');
        setNumSelected(0)
        scatterplot.set({
            opacityBy: 'density',
            cameraView: cview,
        })
        scatterplot.refresh()
    }

    return (
        <div className="select-none fixed bottom-0 rounded-lg right-0 w-fit bg-gray-100/80 p-1 flex text-sm text-gray-500"
        >
            {numSelected > 0 &&
                <p>{numSelected} selected</p>
            }
            <abbr title="hold shift to draw a selection, Shift+Ctrl adds to the current selection, double click on the background removes the selection">
                <InformationCircleIcon className="w-5 h-5 mx-1" />
            </abbr>
            <p>{numPoints} points, {datasetInfo}</p>
        </div>
    )

}