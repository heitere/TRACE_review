/** @type {import('next').NextConfig} */
// const nextConfig = {}
// module.exports = nextConfig


module.exports = () => {
    const rewrites = () => {
      return [
        {
          source: "/backend/:path*",
          destination: "http://backend:8000/backend/:path*",

          // use localhost when running without docker
          //destination: "http://localhost:8000/backend/:path*",
        },
      ];
    };
    return {
      rewrites,
    };
  };

