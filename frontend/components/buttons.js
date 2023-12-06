import { ArrowsPointingInIcon } from "@heroicons/react/24/solid"
import { AdjustmentsHorizontalIcon, MagnifyingGlassIcon } from "@heroicons/react/24/solid"
import { useState, useLayoutEffect, useRef} from "react";

export function MyButton(props) {

  const { onClick, children } = props;

  return (
    <button
      type="button"
      className="h-fit select-none inline-flex justify-center	items-center rounded-md bg-white mx-1 my-1 px-5 pb-2 pt-2.5 text-sm font-medium leading-normal text-gray-700 border-solid drop-shadow-sm border border-gray-600 hover:no-underline hover:bg-gray-100 hover:opacity-75 focus:opacity-100 focus:shadow-none focus:outline-none"
      onClick={onClick}
    >
      {children}
    </button>
  )
}


export function AsyncButton(props) {

  const { onClick, children } = props;
  const [isLoading, setLoading] = useState(false);
  const [buttonWidth, setWidth] = useState("0");
  const ref = useRef();

  const test_dimensions = () => {
    if (ref.current) {
      setWidth((ref.current.offsetWidth + 1).toString());
    }
  }

  useLayoutEffect(() => {
    test_dimensions();
  }, []);

  const handleClick = () => {
    setLoading(true);
    onClick()
      .then(res => setLoading(false));
  }

  return (
    <>
    <button
      type="button"
      className={`min-w-[130px] h-fit select-none justify-center items-center rounded-md bg-white mx-1 my-1 px-5 pb-2 pt-2.5 text-sm font-medium leading-normal text-gray-700 border-solid drop-shadow-sm border border-gray-600 hover:no-underline hover:bg-gray-100 hover:opacity-75 focus:opacity-100 focus:shadow-none focus:outline-none disabled:opacity-75 disabled:bg-gray-100`}
      onClick={handleClick}
      disabled={isLoading}
      ref={ref}
    >
      {isLoading ? "loading..." : children}
    </button>
    </>
  )

}

export function ChevronButton(props) {
  const { innerText, onChange, ...other } = props;
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
    <>
      <MyButton onClick={toggleVisibility}>
        {<>
          <p className="select-none">{innerText}</p>
          {visibility == "visible" ? (<ChevronUpIcon className="w-4 h-4 ml-2" />) : (<ChevronDownIcon className="w-4 h-4 ml-2" />)}
        </>
        }
      </MyButton>
    </>
  )
}

export function ResetButton(props) {
  return (
    <button
      type="button"
      className="select-none inline-block rounded-full bg-gray-600 mx-2 p-2 text-sm font-medium leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 focus:shadow-none focus:outline-none"
      onClick={props.onClick}
    >
      <ArrowsPointingInIcon className="h-6 w-6" />
    </button>
  )
}

export function ZoomButton(props) {
  return (
    <button
      type="button"
      className="select-none inline-block rounded-full bg-gray-600 mx-2 p-2 text-sm font-medium leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 focus:shadow-none focus:outline-none"
      onClick={props.onClick}
    >
      <MagnifyingGlassIcon className="h-6 w-6" />
    </button>
  )
}

export function SettingsButton(props) {
  return (
    <button
      type="button"
      className="select-none inline-block rounded-full bg-gray-600 mx-2 p-2 text-sm font-medium leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 focus:shadow-none focus:outline-none"
      onClick={props.onClick}
    >
      <AdjustmentsHorizontalIcon className="h-6 w-6" />
    </button>
  )
}


export function ChevronRightButton(props) {
  return (
    <button
      type="button"
      className="select-none inline-block rounded-full bg-gray-600 p-2 text-sm font-medium leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 focus:shadow-none focus:outline-none"
      onClick={props.onClick}
    >
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6">
        <path fillRule="evenodd" d="M4.72 3.97a.75.75 0 011.06 0l7.5 7.5a.75.75 0 010 1.06l-7.5 7.5a.75.75 0 01-1.06-1.06L11.69 12 4.72 5.03a.75.75 0 010-1.06zm6 0a.75.75 0 011.06 0l7.5 7.5a.75.75 0 010 1.06l-7.5 7.5a.75.75 0 11-1.06-1.06L17.69 12l-6.97-6.97a.75.75 0 010-1.06z" clipRule="evenodd" />
      </svg>
    </button>
  )
}