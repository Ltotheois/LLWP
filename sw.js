// version 1.1
// change version to reinstall service worker and update Files

var cacheName = 'cache_llwp';
var filesToCache = [
	'/LLWP/',
	'/LLWP/index.html',
	'/LLWP/docs.html',
	'/LLWP/contact.html',
	'/LLWP/download.html',
	'/LLWP/videos.html',
	'/LLWP/LLWP.svg',
	'/LLWP/manifest.json',
	'/LLWP/resources/style.css',
	'/LLWP/resources/download/LLWPvideo.jpg',
	'/LLWP/resources/download/InstallScripts/ubuntu_install.sh',
	'/LLWP/resources/download/InstallScripts/windows_install.bat',
	'/LLWP/resources/fonts/nunito-sans-v6-latin-900.woff2',
	'/LLWP/resources/fonts/nunito-sans-v6-latin-800.woff2',
	'/LLWP/resources/fonts/nunito-sans-v6-latin-700.woff2',
	'/LLWP/resources/fonts/nunito-sans-v6-latin-600.woff2',
	'/LLWP/resources/fonts/nunito-sans-v6-latin-300.woff2',
	'/LLWP/resources/fonts/nunito-sans-v6-latin-700italic.woff2',
	'/LLWP/resources/fonts/nunito-sans-v6-latin-600italic.woff2',
	'/LLWP/resources/fonts/nunito-sans-v6-latin-300italic.woff2'
];

/* Start the service worker and cache all of the app's content */
self.addEventListener('install', function(e) {
  caches.delete(cacheName);
  e.waitUntil(
    caches.open(cacheName).then(function(cache) {
      return cache.addAll(filesToCache);
    })
  );
  self.skipWaiting();
});

/* Serve cached content when offline */
self.addEventListener('fetch', function(e) {
  e.respondWith(
    caches.match(e.request).then(function(response) {
      return response || fetch(e.request);
    })
  );
});