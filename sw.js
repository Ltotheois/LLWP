var cacheName = 'cache_llwp';
var filesToCache = [
	'/LLWP/',
	'/LLWP/index.html',
	'/LLWP/docs.html',
	'/LLWP/contact.html',
	'/LLWP/download.html',
	'/LLWP/videos.html',
	'/LLWP/LLWP.svg',
	'/LLWP/sw.js',
	'/LLWP/manifest.json',
	'/LLWP/resources/style.css',
	'/LLWP/resources/download/LLWPvideo.mp4',
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
self.addEventListener('install', function(event) {
	event.waitUntil(
		caches.open(cacheName).then(function(cache) {
			return cache.addAll(filesToCache);
		})
	);
	self.skipWaiting();
});

/* Serve cached content when offline */
self.addEventListener('fetch', function(event) {
	event.respondWith(
		caches.open(cacheName).then(function(cache){
			cache.match(event.request).then(function(response){
				return response || fetch(event.request);
			});
		});
	);
	
	event.waitUntil(
		caches.open(cacheName).then(function(cache){
			fetch(event.request).then(function(response){
				return cache.put(event.request, response);
			});
		};
	);
});