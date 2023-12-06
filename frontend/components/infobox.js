import { useState } from "react";
import { HoverNote } from "./utils";
import { Switch } from '@headlessui/react'


function HoverText({ pointId, selectedPointColor, pointColors }) {
    return (
        <>
            <p>point: {pointId}</p>
            {selectedPointColor != "none" ?
                <p>
                    {selectedPointColor}: {
                        pointColors["type"] === "categorical" ?
                            pointColors["values"][pointId]
                            : pointColors["values"][pointId].toFixed(2)
                    }
                </p> :
                <p></p>
            }
        </>
    )

}

export function Infobox(props) {
    const {
        scatterplot,
        selectedPointColor,
        pointColors,
        ...other } = props

    const [hoverNoteEnabled, setHoverNoteEnabled] = useState(false)
    const [pointOverUnsubscriber, setPointOverUnsubscriber] = useState(null)
    const [pointOutUnsubscriber, setPointOutUnsubscriber] = useState(null)
    const [hoverState, setHoverState] = useState({
        pointId: 0,
        visibility: "hidden",
        position: [0, 0]
    });

    const handlePointOver = (pointId) => {
        setHoverState({
            pointId: pointId,
            visibility: "visible",
            position: scatterplot.getScreenPosition(pointId)
        })
    }

    const handlePointOut = () => {
        setHoverState(prevState => ({
            ...prevState,
            visibility: "hidden"
        }));
    }

    const getHoverColor = (pointId) => {
        if (pointColors["type"] === "categorical") {
            return pointColors["colorMap"][pointColors["values"][pointId]];
        } else {
            return ("#e9e9e9");
        }
    }


    const handleHoverNoteEnabled = (enabled) => {
        if (enabled) {
            setPointOverUnsubscriber(scatterplot.subscribe('pointover', handlePointOver));
            setPointOutUnsubscriber(scatterplot.subscribe('pointout', handlePointOut));
        }
        else {
            scatterplot.unsubscribe(pointOverUnsubscriber)
            scatterplot.unsubscribe(pointOutUnsubscriber)
        }
        setHoverNoteEnabled(enabled);
    }

    return (
        <>
            <HoverNote
                visible={hoverState.visibility}
                position={hoverState.position}
                color={getHoverColor(hoverState.pointId)}>
                <HoverText pointId={hoverState.pointId} pointColors={pointColors} selectedPointColor={selectedPointColor} />
            </HoverNote>
            <Switch.Group>
                <div className="flex flex-wrap items-center p-3 justify-between">
                    <Switch.Label className="text-sm text-gray-500 w-fit min-w-fit mr-2" htmlFor='hoverSwitch'>point hover info</Switch.Label>
                    <div className="w-1/2 flex items-start justify-start">
                        <Switch
                            id='hoverSwitch'
                            checked={hoverNoteEnabled}
                            onChange={handleHoverNoteEnabled}
                            className={`${hoverNoteEnabled ? 'bg-blue-600' : 'bg-gray-200'
                                } relative inline-flex h-6 w-11 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:ring-offset-2`}
                        >
                            <span
                                className={`${hoverNoteEnabled ? 'translate-x-6' : 'translate-x-1'
                                    } inline-block h-4 w-4 transform rounded-full bg-white transition-transform`}
                            />
                        </Switch>
                    </div>
                </div>
            </Switch.Group>
        </>
    )
}