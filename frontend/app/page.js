import CanvasWrapper from '@/components/canvas_wrapper'
import Scatterplot from '@/components/scatterplot'

export default function Home() {

  return (
    <main className="flex h-screen w-screen flex-row items-top">
      {/* <div className="text-left w-full">
        <h2 className={`m-0 text-2xl font-semibold`}>
          t-SNE Embeddings of human immune cells
        </h2>
      </div> */}
      {/* <div className="mb-0 lg:mb-5 w-full lg:flex text-left">
        <p className={`m-0 text-sm opacity-50`}>
          The dataset from<a className="pointer-events-none lg:pointer-events-auto"
            target="_blank"
            rel="noopener noreferrer"
            href="https://www.nature.com/articles/s41592-021-01336-8"
          >Luecken et al. (2022)</a> contains
          33506 cells from five donors across four sequencing technologies and two tissues (peripheral blood and bone marrow).
        </p>
      </div> */}
      {/* <CanvasWrapper/> */}
      <Scatterplot />
    </main>
  )
}
