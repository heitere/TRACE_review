import { Component, useState } from 'react';
import Select from 'react-select'


class NamedSlider extends Component {
  constructor(props) {
    super(props);
    this.state = {
      value: props.defaultValue,
    };

    this.showCurrentValue = (value) => {
      this.setState({ value: value });
      this.props.onChange(value)
    };
  }

  render() {
    return (
      <>
        <label
          className="select-none text-sm text-gray-500 min-w-fit mr-2">
          {
            this.props.names !== undefined ?
              this.props.label + " " + this.props.names[this.state.value] :
              this.props.label + " " + this.state.value
          }
        </label>
        <input
          className="transparent h-[3px] w-1/2 min-w-[100px] cursor-pointer appearance-none border-transparent bg-neutral-300 my-2"
          type="range" min={this.props.min} max={this.props.max} step={this.props.step} defaultValue={this.props.defaultValue}
          onChange={(event) => this.showCurrentValue(+event.target.value)} />
      </>
    );
  }
}

export default NamedSlider



const RadioButton = ({ onChange, value, checked }) => (
  <label>
    <input type="radio" name="radio-button-group" value={value} onChange={onChange} checked={checked} /> {value}
  </label>
);


export function RadioSelect(props) {
  const { onChange, defaultValue, options, ...other } = props;
  const [currentValue, setCurrentValue] = useState(defaultValue);

  function onRadioChange(event) {
    setCurrentValue(event.target.value);
    onChange(event.target.value);
  }

  return (
    <div className='flex flex-col w-fit h-fit'>
      {
        options.map(option =>
          <RadioButton key={option} value={option} checked={option == currentValue} onChange={onRadioChange} />
        )
      }
    </div>
  )
}


export function DropdownSelectOld(props) {
  const { options, selected, onChange, id, ...other } = props

  return (
    <select id={id} value={selected} onChange={(event) => onChange(event.target.value)}>
      {options.map((option, optionIdx) => (
        <option key={optionIdx} value={option}>
          {option}
        </option>
      ))}
    </select>
  );
}


export function ReactSelect(props) {
  const { options, selected, onChange, menuPlacement = 'auto', ...other } = props

  return (
    <Select
      className="w-1/2 min-w-max text-slate-600 text-left"
      options={options}
      isClearable={false}
      isSearchable={true}
      menuPlacement={menuPlacement}
      value={{ 'value': selected, 'label': selected }}
      onChange={(selection) => { onChange(selection["value"]) }}
    />
  )
}



export function HoverNote(props) {
  const { visible, color, position, children } = props;
  return (
    <button
      type="button"
      className={"select-none fixed rounded-md bg-white/90 outline outline-2 text-left m-2 p-2 text-base font-medium leading-normal text-black"}
      style={{
        'outlineColor': color,
        'visibility': visible,
        'top': position[1],
        'left': position[0],
      }}
    >
      {children}
    </button >
  )
}